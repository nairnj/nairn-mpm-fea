/********************************************************************************
    Viscoelastic.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Feb 5 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/Viscoelastic.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "System/UnitsController.hpp"
#include "Exceptions/MPMWarnings.hpp"

// class statics
int Viscoelastic::warnExcessiveX = -1;

#pragma mark Viscoelastic::Constructors and Destructors

// Constructors
Viscoelastic::Viscoelastic()
{
}

// Constructors
Viscoelastic::Viscoelastic(char *matName) : MaterialBase(matName)
{
	ntaus=-1;
    Gk=NULL;
    tauk=NULL;
    G0=0.;
    currentGk=0;
    currentTauk=0;
	aI=40.;
	
	// MGEOS variables
	pressureLaw = LINEAR_PRESSURE;
	mptrHistory=-1;		// Store J and J res if using MGEOS
	gamma0=1.64;		// dimensionless
	C0=4004000.;		// m/sec
	S1=1.35;			// dimsionless
	S2=0.;				// dimsionless
	S3=0.;				// dimsionless
	Kmax=20.;			// maxium relative increase allows in K
}

#pragma mark Viscoelastic::Initialization

// print mechanical properties to output window
void Viscoelastic::PrintMechanicalProperties(void) const
{
	if(pressureLaw==LINEAR_PRESSURE)
    {	cout << "Pressure law: linear elastic" << endl;
		PrintProperty("K",K*UnitsController::Scaling(1.e-6),"");
		PrintProperty("a",aI,"");
		cout << endl;
	}
	else
	{	// core properties
		cout << "Pressure law: MG-EOS" << endl;
		PrintProperty("C0",C0*UnitsController::Scaling(1.e-3),UnitsController::Label(ALTVELOCITY_UNITS));
		PrintProperty("gam0",gamma0,"");
		PrintProperty("K",K*UnitsController::Scaling(1.e-6),"");
		cout << endl;
		
		PrintProperty("S1",S1,"");
		PrintProperty("S2",S2,"");
		PrintProperty("S3",S3,"");
		cout << endl;
		
		// effective volumetric CTE (in ppm/K) alpha = rho0 gamma0 Cv / K
		PrintProperty("aI",1.e6*CTE,"");
		PrintProperty("T0",thermal.reference,"K");
		cout <<  endl;
		
		// Kmax
		PrintProperty("Kmax",Kmax," K0");
		PrintProperty("Xmax",Xmax,"");
		cout <<  endl;
	}
	PrintProperty("G0",G0*UnitsController::Scaling(1.e-6),"");
	PrintProperty("ntaus",(double)ntaus,"");
    cout <<  endl;
	
	int i;
    for(i=0;i<ntaus;i++)
	{	PrintProperty("i",(double)(i+1),"");
		PrintProperty("Gk",Gk[i]*UnitsController::Scaling(1.e-6),"");
		PrintProperty("tauk",tauk[i],UnitsController::Label(TIME_UNITS));
        cout << endl;
    }
	
	// For information, E and nu at time 0
	PrintProperty("E(0)",9.*Kered*Gered*rho/(3.*Kered+Gered)*UnitsController::Scaling(1.e-6),"");
	PrintProperty("nu(0)",(3.*Kered-2.*Gered)/(6.*Kered+2.*Gered),"");
	cout << endl;
}
    
// Read material properties
char *Viscoelastic::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"G0")==0)
		return UnitsController::ScaledPtr((char *)&G0,gScaling,1.e6);
    
    else if(strcmp(xName,"K")==0)
		return UnitsController::ScaledPtr((char *)&K,gScaling,1.e6);
    
    else if(strcmp(xName,"alpha")==0)
        return((char *)&aI);
    
    else if(strcmp(xName,"ntaus")==0)
    {	input=INT_NUM;
        return((char *)&ntaus);
    }
    
    else if(strcmp(xName,"Gk")==0)
    {	if(Gk==NULL)
        {   if(ntaus<=0)
				ThrowSAXException("Gk found before number of taus specified.");
            Gk=new double[ntaus];
        }
        currentGk++;
        if(currentGk>ntaus)
			ThrowSAXException("Too many Gk's given.");
		return UnitsController::ScaledPtr((char *)&Gk[currentGk-1],gScaling,1.e6);
    }
    
    else if(strcmp(xName,"tauk")==0)
    {	if(tauk==NULL)
        {   if(ntaus<=0)
				ThrowSAXException("tauk found before number of taus specified.");
            tauk=new double[ntaus];
        }
        currentTauk++;
        if(currentTauk>ntaus)
			ThrowSAXException("Too many tauk's given.");
        return((char *)&tauk[currentTauk-1]);
    }
	
	// The follow allow use of MG-EOS law for pressure dependence
	
    else if(strcmp(xName,"pressureLaw")==0)
    {	input=INT_NUM;
        return((char *)&pressureLaw);
    }
    
    else if(strcmp(xName,"gamma0")==0)
    {	input=DOUBLE_NUM;
        return((char *)&gamma0);
    }
    
    else if(strcmp(xName,"C0")==0)
    {	input=DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&C0,gScaling,1.e3);
    }
    
    else if(strcmp(xName,"S1")==0)
    {	input=DOUBLE_NUM;
        return((char *)&S1);
    }
	
    else if(strcmp(xName,"S2")==0)
    {	input=DOUBLE_NUM;
        return((char *)&S2);
    }
	
    else if(strcmp(xName,"S3")==0)
    {	input=DOUBLE_NUM;
        return((char *)&S3);
    }

    else if(strcmp(xName,"Kmax")==0)
    {	input=DOUBLE_NUM;
        return((char *)&Kmax);
    }
	    
    return(MaterialBase::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
const char *Viscoelastic::VerifyAndLoadProperties(int np)
{
	// check properties
    if(currentGk<ntaus || currentTauk<ntaus)
		return "Insufficient Gk or tauk for expected number of taus.";
    
    if(ntaus<0)
		return "Number of taus was never entered.";
    
    if(K<0)
		return "Required bulk modulus not given.";
    
    // zero time shear modulus
    Gered = G0;
	TwoGkred = new double[ntaus];
    for(int k=0;k<ntaus;k++)
    {   Gered += Gk[k];
		TwoGkred[k] = 2.*Gk[k]/rho;
    }
	
	// Convert to specific moduli
	Gered /= rho;
	
	// to absolute CTE and CME
	CTE = 1.e-6*aI;
	CME = betaI*concSaturation;
	
	// bulk modulus
	if(pressureLaw==LINEAR_PRESSURE)
	{	Kered = K/rho;
	}
	
	else if(pressureLaw==MGEOS_PRESSURE)
	{	// Use in place of C0^2. Units are L^2/sec^2 = F/L^2 L^3/mass
		// Equal to reduced bulk modulus
		Kered = C0squared = C0*C0;
		
		// Initial bulk modulus
		K = rho*C0squared;
		
		// expansion coefficients - affect on pressure is handled by eos
		// find for printing (and maybe future large deformation shear)
		double effAlpha = (heatCapacity*gamma0)/C0squared;
		CTE = effAlpha/3.;
		
		// this material not coupled to moisture expansion
		betaI = 0.;
		CME = 0.;
		
		// warning
		if(warnExcessiveX<0)
			warnExcessiveX = warnings.CreateWarning("compressive strain has exceeded MGEOS law range in Viscoleastic material",-1,0);
		Xmax = GetMGEOSXmax(gamma0,S1,S2,S3,Kmax);
	}
	
	else
		return "Invalid option for the pressure law";
	
    // for Cp-Cv (units nJ/(g-K^2))
    Ka2sp = Kered*CTE*CTE;
	
	// call super class
	return MaterialBase::VerifyAndLoadProperties(np);
}

// plane stress not allowed in viscoelasticity
void Viscoelastic::ValidateForUse(int np) const
{	if(np==PLANE_STRESS_MPM)
	{	throw CommonException("Viscoelastic materials are not allowed in plane stress calculations",
							  "Viscoelastic::ValidateForUse");
	}
	
	//call super class (why can't call super class?)
	MaterialBase::ValidateForUse(np);
}

#pragma mark Viscoelastic::History Data Methods

// create and return pointer to history variables
// initialize all to zero
char *Viscoelastic::InitHistoryData(char *pchr,MPMBase *mptr)
{
	// if none, only need particle history
#ifdef USE_KIRCHOFF_STRESS
    if(ntaus==0)
	{	char *p;
		if(pchr==NULL)
			p = new char[sizeof(double *)];
		else
			p = pchr;
		double **h = (double **)p;
		mptrHistory=0;
		h[mptrHistory] = new double[2];
		h[mptrHistory][MGJ_HISTORY] = 1.;					// J
		h[mptrHistory][MGJRES_HISTORY] = 1.;				// Jres (MGEOS only)
		return p;
	}
#else
    if(ntaus==0)
	{	if(pressureLaw==MGEOS_PRESSURE)
		{	char *p;
			if(pchr==NULL)
				p = new char[sizeof(double *)];
			else
				p = pchr;
			double **h = (double **)p;
			mptrHistory=0;
			h[mptrHistory] = new double[2];
			h[mptrHistory][MGJ_HISTORY] = 1.;					// J
			h[mptrHistory][MGJRES_HISTORY] = 1.;				// Jres
			return p;
		}
		else
			return NULL;
	}
#endif
    
    // allocate array of double pointers (3)
	int blocks;
	if(fmobj->IsThreeD())
		blocks = 6;
	else
		blocks = 4;
	
	// one extra for any additional history variables
#ifdef USE_KIRCHOFF_STRESS
	mptrHistory=blocks;
	blocks++;
#else
	if(pressureLaw==MGEOS_PRESSURE)
	{	mptrHistory=blocks;
		blocks++;
	}
#endif
	
	// history variables are pointers to arrays of doubles
    char *p = new char[sizeof(double *)*blocks];
	
    double **h = (double **)p;
    
    // for each allocate ntaus doubles
    //	h[ij_HISTORY][0] to h[ij_HISTORY][ntaus-1] can be read
    //	in MPMConstitutiveLaw() by casting mptr->GetHistoryPtr(0) pointer as
    //		double **h=(double **)(mptr->GetHistoryPtr(0))
    h[XX_HISTORY] = new double[ntaus];
    h[YY_HISTORY] = new double[ntaus];
    h[XY_HISTORY] = new double[ntaus];
	h[ZZ_HISTORY] = new double[ntaus];
	if(fmobj->IsThreeD())
	{	h[XZ_HISTORY] = new double[ntaus];
		h[YZ_HISTORY] = new double[ntaus];
	}		
    
    // initialize to zero
    int k;
    for(k=0;k<ntaus;k++)
    {	h[XX_HISTORY][k] = 0.;
        h[YY_HISTORY][k] = 0.;
        h[XY_HISTORY][k] = 0.;
		h[ZZ_HISTORY][k] = 0.;
		if(fmobj->IsThreeD())
		{	h[XZ_HISTORY][k] = 0.;
			h[YZ_HISTORY][k] = 0.;
		}
    }
    
	// extra particle history variables
#ifdef USE_KIRCHOFF_STRESS
	h[mptrHistory] = new double[2];
	h[mptrHistory][MGJ_HISTORY] = 1.;					// J
	h[mptrHistory][MGJRES_HISTORY] = 1.;				// Jres (MGEOS only)
#else
	if(mptrHistory>0)
	{	h[mptrHistory] = new double[2];
		h[mptrHistory][MGJ_HISTORY] = 1.;					// J
		h[mptrHistory][MGJRES_HISTORY] = 1.;				// Jres
	}
#endif
	
    return p;
}

// Get J and Jres only and only in nonlinear pressure law
double Viscoelastic::GetHistory(int num,char *historyPtr) const
{
    double history=0.;
#ifdef USE_KIRCHOFF_STRESS
	if(num>0 && num<=MGJRES_HISTORY+1)
	{	double **h =(double **)historyPtr;
		history = h[mptrHistory][num-1];
	}
#else
	if(num>0 && num<=MGJRES_HISTORY+1)
	{	if(pressureLaw==LINEAR_PRESSURE)
			history = 1.;
		else
		{	double **h =(double **)historyPtr;
			history = h[mptrHistory][num-1];
		}
	}
#endif
    return history;
}

#pragma mark Viscoelastic::Methods

/* Take increments in strain and calculate new Particle: strains, rotation strain,
	stresses, strain energy,
	du are (gradient rates X time increment) to give deformation gradient change
	For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMETRIC_MPM, otherwise dvzz=0
	This material tracks pressure and stores deviatoric stress only in particle stress tensor
 */
void Viscoelastic::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res,int historyOffset) const
{
	// current previous deformation gradient and stretch
	Matrix3 pFnm1 = mptr->GetDeformationGradientMatrix();
	
    // get incremental deformation gradient and decompose it
	const Matrix3 dF = du.Exponential(incrementalDefGradTerms);
    Matrix3 dR;
    Matrix3 dVstretch = dF.LeftDecompose(&dR,NULL);
	
	// decompose to get previous stretch
	Matrix3 Vnm1 = pFnm1.LeftDecompose(NULL,NULL);
	
	// get strain increments de = (dV-I) dR Vnma dRT
	dVstretch(0,0) -= 1.;
	dVstretch(1,1) -= 1.;
	dVstretch(2,2) -= 1.;
	Matrix3 Vrot = Vnm1.RMRT(dR);
	Matrix3 detot = dVstretch*Vrot;
	
	// Update total deformation gradient
	Matrix3 pF = dF*pFnm1;
	mptr->SetDeformationGradientMatrix(pF);
	
    // Effective strain by deducting thermal strain (no shear thermal strain because isotropic)
	double eres=CTE*res->dT;
	if(DiffusionTask::active)
		eres+=CME*res->dC;
	
	// update pressure
	double dTq0 = 0.,dispEnergy = 0.;
	double traceDe = detot.trace();
	double delV = traceDe - 3.*eres;
#ifdef USE_KIRCHOFF_STRESS
	// tracking J
	double detdF = dF.determinant();
	double **h =(double **)(mptr->GetHistoryPtr(0));
	double J = detdF*h[mptrHistory][MGJ_HISTORY];
	h[mptrHistory][MGJ_HISTORY] = J;
	UpdatePressure(mptr,delV,res,eres,detdF,J,delTime,dTq0,dispEnergy);
#else
	UpdatePressure(mptr,delV,res,eres,&dF,delTime,dTq0,dispEnergy);
#endif
	
	// deviatoric strains increment in de
	// Actually de is finding 2*(dev e) to avoid many multiplies by two
	Tensor de;
	double dV = traceDe/3.;
	de.xx = 2.*(detot(0,0) - dV);
	de.yy = 2.*(detot(1,1) - dV);
	de.zz = 2.*(detot(2,2) - dV);
	de.xy = 2.*detot(0,1);
	if(np==THREED_MPM)
	{	de.xz = 2.*detot(0,2);
		de.yz = 2.*detot(1,2);
	}
	
	// Find initial 2*e(t) (deviatoric strain) in ed
	Tensor ed;
	double thirdV = Vrot.trace()/3.;
	ed.xx = 2.*(Vrot(0,0)-thirdV);
	ed.yy = 2.*(Vrot(1,1)-thirdV);
	ed.zz = 2.*(Vrot(2,2)-thirdV);
	ed.xy = 2.*Vrot(0,1);
	if(np==THREED_MPM)
	{	ed.xz = 2.*Vrot(0,2);
		ed.yz = 2.*Vrot(1,2);
	}
	
	// increment particle deviatoric stresses - elastic part
	double dsig[6];
	dsig[XX] = Gered*de.xx;
	dsig[YY] = Gered*de.yy;
	dsig[ZZ] = Gered*de.zz;
	dsig[XY] = Gered*de.xy;
	if(np==THREED_MPM)
	{	dsig[XZ] = Gered*de.xz;
		dsig[YZ] = Gered*de.yz;
	}
	
	// get internal variable increments, update them, add to incremental stress, and get dissipated energy6
	Tensor dak;
    double **ak =(double **)(mptr->GetHistoryPtr(0));
	int k;
    for(k=0;k<ntaus;k++)
    {   double tmp = exp(-delTime/tauk[k]);
		double tmpm1 = tmp-1.;
		double tmpp1 = tmp+1.;
		double arg = 0.25*delTime/tauk[k];					// 0.25 because e's have factor of 2
		dak.xx = tmpm1*ak[XX_HISTORY][k] + arg*(tmpp1*ed.xx + de.xx);
		dak.yy = tmpm1*ak[YY_HISTORY][k] + arg*(tmpp1*ed.yy + de.yy);
		dak.xy = tmpm1*ak[XY_HISTORY][k] + arg*(tmpp1*ed.xy + de.xy);
		dak.zz = tmpm1*ak[ZZ_HISTORY][k] + arg*(tmpp1*ed.zz + de.zz);
		ak[XX_HISTORY][k] += dak.xx;
		ak[YY_HISTORY][k] += dak.yy;
		ak[ZZ_HISTORY][k] += dak.zz;
		ak[XY_HISTORY][k] += dak.xy;
		
		// add to stress increments
		dsig[XX] -= TwoGkred[k]*dak.xx;
		dsig[YY] -= TwoGkred[k]*dak.yy;
		dsig[ZZ] -= TwoGkred[k]*dak.zz;
		dsig[XY] -= TwoGkred[k]*dak.xy;
		
		// dissipation
		dispEnergy += TwoGkred[k]*(dak.xx*(0.5*(ed.xx+0.5*de.xx)-ak[XX_HISTORY][k]+0.5*dak.xx)
								   + dak.yy*(0.5*(ed.yy+0.5*de.yy)-ak[YY_HISTORY][k]+0.5*dak.yy)
								   + dak.zz*(0.5*(ed.zz+0.5*de.zz)-ak[ZZ_HISTORY][k]+0.5*dak.zz)
								   + dak.xy*(0.5*(ed.xy+0.5*de.xy)-ak[XY_HISTORY][k]+0.5*dak.xy));
		
		// extra terms for 3D
		if(np==THREED_MPM)
		{	// internal variables
			dak.xz = tmpm1*ak[XZ_HISTORY][k] + arg*(tmpp1*ed.xz + de.xz);
			dak.yz = tmpm1*ak[YZ_HISTORY][k] + arg*(tmpp1*ed.yz + de.yz);
			ak[XZ_HISTORY][k] += dak.xz;
			ak[YZ_HISTORY][k] += dak.yz;
			
			// stresses
			dsig[XZ] -= TwoGkred[k]*dak.xz;
			dsig[YZ] -= TwoGkred[k]*dak.yz;
			
			// dissipation
			dispEnergy += TwoGkred[k]*(dak.xz*(0.5*(ed.xz+0.5*de.xz)-ak[XZ_HISTORY][k]+0.5*dak.xz)
									   + dak.yz*(0.5*(ed.yz+0.5*de.yz)-ak[YZ_HISTORY][k]+0.5*dak.yz));
		}
	}
	
	// Update particle deviatoric stresses
	Tensor *sp=mptr->GetStressTensor();
	//Tensor st0 = *sp;
	if(np==THREED_MPM)
	{	// incremental rotate of prior stress
		Matrix3 stn(sp->xx,sp->xy,sp->xz,sp->xy,sp->yy,sp->yz,sp->xz,sp->yz,sp->zz);
		Matrix3 str = stn.RMRT(dR);

#ifdef USE_KIRCHOFF_STRESS
		// convert sigma(n)/rho(n) to sigma(n)/rho(n+1) and add dsigma/rho(n+1)
		sp->xx = detdF*str(0,0)+J*dsig[XX];
		sp->yy = detdF*str(1,1)+J*dsig[YY];
		sp->xy = detdF*str(0,1)+J*dsig[XY];
		sp->zz = detdF*str(2,2)+J*dsig[ZZ];
		sp->yz = detdF*str(1,2)+J*dsig[YZ];
		sp->xz = detdF*str(0,2)+J*dsig[XZ];
#else
		sp->xx = str(0,0)+dsig[XX];
		sp->yy = str(1,1)+dsig[YY];
		sp->xy = str(0,1)+dsig[XY];
		sp->zz = str(2,2)+dsig[ZZ];
		sp->yz = str(1,2)+dsig[YZ];
		sp->xz = str(0,2)+dsig[XZ];
#endif
	}
	else
	{	// incremental rotate of prior stress
		Matrix3 stn(sp->xx,sp->xy,sp->xy,sp->yy,sp->zz);
		Matrix3 str = stn.RMRT(dR);
		
#ifdef USE_KIRCHOFF_STRESS
		// convert sigma(n)/rho(n) to sigma(n)/rho(n+1) and add dsigma/rho(n+1)
		sp->xx = detdF*str(0,0)+J*dsig[XX];
		sp->yy = detdF*str(1,1)+J*dsig[YY];
		sp->xy = detdF*str(0,1)+J*dsig[XY];
		sp->zz = detdF*sp->zz+J*dsig[ZZ];
#else
		sp->xx = str(0,0)+dsig[XX];
		sp->yy = str(1,1)+dsig[YY];
		sp->xy = str(0,1)+dsig[XY];
		sp->zz += dsig[ZZ];
#endif
	}
	
	// incremental work energy = shear energy (dilation and residual energy done in update pressure)
    double shearEnergy = sp->xx*detot(0,0) + sp->yy*detot(1,1) + sp->zz*detot(2,2) + sp->xy*de.xy;
    if(np==THREED_MPM)
    {   shearEnergy += sp->xz*de.xz + sp->yz*de.yz;
    }
    mptr->AddWorkEnergyAndResidualEnergy(shearEnergy,0.);
	
    // disispated energy per unit mass (dPhi/(rho0 V0)) (nJ/g)
    mptr->AddPlastEnergy(dispEnergy);
    
    // heat energy is Cv dT - dPhi
    // The dPhi is subtracted here because it will show up in next
    //		time step within Cv dT (if adibatic heating occurs)
    // The Cv dT was done in update pressure
    IncrementHeatEnergy(mptr,res->dT,dTq0,dispEnergy);
}

// This method handles the pressure equation of state. Its tasks are
// 1. Calculate the new pressure
// 2. Update particle pressure
// 3. Increment the particle energy
#ifdef USE_KIRCHOFF_STRESS
void Viscoelastic::UpdatePressure(MPMBase *mptr,double delV,ResidualStrains *res,double eres,double detdF,double J,
								  double delTime,double &dTq0,double &AVEnergy) const
#else
void Viscoelastic::UpdatePressure(MPMBase *mptr,double delV,ResidualStrains *res,double eres,const Matrix3 *dF,
								  double delTime,double &dTq0,double &AVEnergy) const
#endif
{
	if(pressureLaw==LINEAR_PRESSURE)
	{	// pressure change
#ifdef USE_KIRCHOFF_STRESS
		double dP = (detdF-1.)*mptr->GetPressure()-J*Kered*delV;
#else
		double dP = -Kered*delV;
#endif
		
		// artifical viscosity
		// delV is total incremental volumetric strain = total Delta(V)/V
		double QAVred = 0.;
		if(delV<0. && artificialViscosity)
		{	// Wants K/rho
#ifdef USE_KIRCHOFF_STRESS
			QAVred = GetArtificalViscosity(delV/delTime,sqrt(Kered*J),mptr);
#else
			QAVred = GetArtificalViscosity(delV/delTime,sqrt(Kered),mptr);
#endif
			if(ConductionTask::AVHeating) AVEnergy += fabs(QAVred*delV);
		}
		
		// increment pressure
		mptr->IncrementPressure(dP+QAVred);
		
		// work energy is dU = -P dV + s.de(total)
		// Here do hydrostatic term
		// Work energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
		double avgP = mptr->GetPressure()-0.5*dP;
		mptr->AddWorkEnergyAndResidualEnergy(-avgP*delV,-3.*avgP*eres);
		
		// Isoentropic temperature rise = -(K alpha T)/(rho Cv) (dV/V)
		dTq0 += -Kered*CTE*mptr->pPreviousTemperature*(delV+3*eres)/GetHeatCapacity(mptr);;
	}
	else
	{	// J is total volume change - may need to reference to free-swelling volume if that works
		// Note that swelling looks like a problem because the sums of strains needs to be adjusted
		//		to stress-free state
		
		// history pointer
		double **h =(double **)(mptr->GetHistoryPtr(0));
		
#ifndef USE_KIRCHOFF_STRESS
		// tracking J
		double detdF = dF->determinant();
		double J = detdF*h[mptrHistory][MGJ_HISTORY];
		h[mptrHistory][MGJ_HISTORY] = J;
#endif

		// account for residual strains (needed in tension, but must always track)
		double dJres = exp(3.*CTE*res->dT);
		double Jres = dJres*h[mptrHistory][MGJRES_HISTORY];
		h[mptrHistory][MGJRES_HISTORY] = Jres;
		
		// previous pressure
		double P,P0 = mptr->GetPressure();
		double delVMG = 1. - 1./detdF;			// (Vnp1-Vn)/Vnp1
		double Kred;
		
		// M-G EOS
		// Want specific pressure or pressure over current density (using J = rho0/rho)
		if(J<1.)
		{	// new compression J(k+1) = 1-x(k+1)
			double x = 1.-J;
			
			if(x<=Xmax)
			{	// compression law
				// denominator = 1 - S1*x - S2*x^2 - S3*x^3
				double denom = 1./(1. - x*(S1 + x*(S2 + x*S3)));
				
				// current effective and reduced (by rho0) bulk modulus
				Kred = C0squared*(1.-0.5*gamma0*x)*denom*denom;
			}
			else
			{	// law not valid if gets too high
				if(warnings.Issue(warnExcessiveX,-1)==GAVE_WARNING)
				{
#pragma omp critical (output)
					{	cout << "# Excessive x = " << x << " causing Kred to increase more than " << Kmax << " fold" << endl;
						mptr -> Describe();
					}
				}
				
				// truncate effective and reduced (by rho0) bulk modulus
				Kred = C0squared*Kmax;
			}
			
			// Pressure from bulk modulus and an energy term
			double e = mptr->GetInternalEnergy();
#ifdef USE_KIRCHOFF_STRESS
			P = J*(Kred*x + gamma0*e);
#else
			P = Kred*x + gamma0*e;
#endif
			
			// particle isentropic temperature increment
			dTq0 += -J*gamma0*mptr->pPreviousTemperature*delVMG;
		}
		else
		{   // In tension hyperelastic law P = - K0(J-1)
			double Jeff = J/Jres;
			Kred = C0squared*Jeff;
#ifdef USE_KIRCHOFF_STRESS
			P = -J*C0squared*(Jeff-1.);
#else
			P = -C0squared*(Jeff-1.);
#endif
			
			// particle isentropic temperature increment
			double Kratio = Jeff;
			dTq0 += -J*Kratio*gamma0*mptr->pPreviousTemperature*delVMG;
		}
		
		// artifical viscosity
		// delVMG is total incremental volumetric strain = total Delta(V)/V
		double QAVred = 0.;
		if(delVMG<0. && artificialViscosity)
		{	QAVred = GetArtificalViscosity(delVMG/delTime,sqrt(Kred*J),mptr);
			if(ConductionTask::AVHeating) AVEnergy += fabs(QAVred*delVMG);
		}
		
		// set final pressure
		mptr->SetPressure(P+QAVred);
		
		// work energy is dU = -P dV + s.de(total)
		// Here do hydrostatic terms, deviatoric later
		double avgP = 0.5*(P0+P);
		double delVres = 1. - 1./dJres;
		mptr->AddWorkEnergyAndResidualEnergy(-avgP*delVMG,-avgP*delVres);
		
	}
}

// convert J to K using isotropic method
Vector Viscoelastic::ConvertJToK(Vector d,Vector C,Vector J0,int np)
{	double nuLS = (3.*Kered-2.*Gered)/(6.*Kered+2.*Gered);
	return IsotropicJToK(d,C,J0,np,nuLS,Gered*rho);
}

// From thermodyanamics Cp-Cv = 9 K a^2 T/rho
// Ka2sp in nJ/(g-K^2) so output in nJ/(g-K)
// Here using K0 and rho0 - could modify if needed
double Viscoelastic::GetCpMinusCv(MPMBase *mptr) const
{   return mptr!=NULL ? Ka2sp*mptr->pPreviousTemperature : Ka2sp*thermal.reference;
}

#pragma mark Viscoelastic::Accessors

// return material type
const char *Viscoelastic::MaterialType(void) const { return "Viscoelastic"; }

// Return the material tag
int Viscoelastic::MaterialTag(void) const { return VISCOELASTIC; }

// Calculate wave speed in mm/sec
// Uses sqrt((K +4Ge/3)/rho) which is probably the maximum wave speed possible
double Viscoelastic::WaveSpeed(bool threeD,MPMBase *mptr) const { return sqrt((Kered + 4.*Gered/3.)); }

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatory stress into full stress
Tensor Viscoelastic::GetStress(Tensor *sp,double pressure,MPMBase *mptr) const
{	Tensor stress = *sp;
    stress.xx -= pressure;
    stress.yy -= pressure;
    stress.zz -= pressure;
    return stress;
}

// if a subclass material supports artificial viscosity, override this and return TRUE
bool Viscoelastic::SupportsArtificialViscosity(void) const { return true; }

// Calculate current wave speed in mm/sec. Uses sqrt((K+4G/3)/rho) which is dilational wave speed
// but K/rho = Kred*J and G/rho = Gred*J (in mm^2/sec^2)
double Viscoelastic::CurrentWaveSpeed(bool threeD,MPMBase *mptr,int offset) const
{
	double KcurrRed = Kered;
	
	if(pressureLaw==MGEOS_PRESSURE)
    {	// compressive volumetric strain x = 1-J
		double J = mptr->GetRelativeVolume();
		
		// get K/rho0, but this ignores slope of energy term
		double KcurrRed;
		if(J<1.)
		{   double x = 1. - J;
			
			if(x<Xmax)
			{	// compression law
				// denominator = 1 - S1*x - S2*x^2 - S3*x^3
				double denom = 1./(1. - x*(S1 + x*(S2 + x*S3)));
			
				// current effective and reduced (by rho0) bulk modulus
				KcurrRed = C0squared*(1.-0.5*gamma0*x)*denom*denom;
			}
			else
			{	// truncate if law seems bad
				KcurrRed = C0squared*Kmax;
			}
		}
		
		//KcurrRed *= J;          // converts to K/rho, but keep K/rho0 for hypoelastic material
	}
    
    // return current save speed
    return sqrt((KcurrRed + 4.*Gered/3.));
}

#ifdef USE_KIRCHOFF_STRESS
// Get current relative volume change = J (which this material tracks)
double Viscoelastic::GetCurrentRelativeVolume(MPMBase *mptr,int offset) const
{   double **h =(double **)mptr->GetHistoryPtr(offset);
	return h[mptrHistory][MGJ_HISTORY];
}
#endif

