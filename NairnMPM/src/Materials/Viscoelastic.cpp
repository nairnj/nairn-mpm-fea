/********************************************************************************
    Viscoelastic.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Feb 5 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
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

// Constructor
Viscoelastic::Viscoelastic(char *matName,int matID) : MaterialBase(matName,matID)
{
	ntaus=-1;
    Gk=NULL;
    tauk=NULL;
    G0=0.;
    currentGk=0;
    currentTauk=0;
	aI=40.;
    
    ntausK=0;
    Kk=NULL;
    tauKk=NULL;
    currentKk=0;
    currentTauKk=0;
	
	// MGEOS variables
	pressureLaw = LINEAR_PRESSURE;
	gamma0=1.64;		// dimensionless
	C0=4004000.;		// m/sec
	S1=1.35;			// dimsionless
	S2=0.;				// dimsionless
	S3=0.;				// dimsionless
	Kmax=-1.;			// maxium relative increase allows in K (default no limit)
    
    // WLF parameters
    Tref = -1.;
    C1base10 = 17.44;
    C2 = 51.6;
    
    // Moisture parameters
    mref = -1.;
    Cm1base10 = 10.;
    Cm2base10 = 0.025/0.4;
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
    else if(pressureLaw==TIME_DEPENDENT_PRESSURE)
    {   cout << "Pressure law: time-dependent bulk modulus" << endl;
        PrintProperty("K(0)",rho*Kered*UnitsController::Scaling(1.e-6),"");
        PrintProperty("K0",K*UnitsController::Scaling(1.e-6),"");
        PrintProperty("ntausK",(double)ntausK,"");
        cout <<  endl;
        
        for(int i=0;i<ntaus;i++)
        {   PrintProperty("  i",(double)(i+1),"");
            PrintProperty("Kk",Kk[i]*UnitsController::Scaling(1.e-6),"");
            PrintProperty("tauKk",tauKk[i],UnitsController::Label(TIME_UNITS));
            cout << endl;
        }
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
		if(Kmax>0.)
		{	PrintProperty("Kmax",Kmax," K0");
			PrintProperty("Xmax",Xmax,"");
		}
		else
			cout << "Kmax= no limit";
		cout <<  endl;
	}
    
    PrintProperty("G(0)",rho*Gered*UnitsController::Scaling(1.e-6),"");
	PrintProperty("G0",G0*UnitsController::Scaling(1.e-6),"");
	PrintProperty("ntaus",(double)ntaus,"");
    cout <<  endl;
	
	int i;
    for(i=0;i<ntaus;i++)
	{	PrintProperty("  i",(double)(i+1),"");
		PrintProperty("Gk",Gk[i]*UnitsController::Scaling(1.e-6),"");
		PrintProperty("tauk",tauk[i],UnitsController::Label(TIME_UNITS));
        cout << endl;
    }
	
	// For information, E and nu at time 0
	PrintProperty("E(0)",9.*Kered*Gered*rho/(3.*Kered+Gered)*UnitsController::Scaling(1.e-6),"");
	PrintProperty("nu(0)",(3.*Kered-2.*Gered)/(6.*Kered+2.*Gered),"");
	cout << endl;
    
    // WLF properties
    if(Tref>=0.)
    {   PrintProperty("Tref",Tref,"K");
        PrintProperty("C1",C1base10,"");
        PrintProperty("C2",C2,"");
        cout << endl;
    }
    else if(mref<0.)
        cout << "Isothermal and isosolvent viscoelasticity" << endl;
    else
        cout << "Isothermal viscoelasticity" << endl;

    // WLF moisture properties
    if(mref>=0.)
    {   PrintProperty("mref",mref*concSaturation,"");
        PrintProperty("Cm1",Cm1base10,"");
        PrintProperty("Cm2",Cm2base10,"");
        cout << endl;
    }
    else if(Tref>=0.)
        cout << "Isosolvent viscoelasticity" << endl;
}
    
// Read material properties
// throws std::bad_alloc, SAXException()
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
				ThrowSAXException("Gk found before number of G taus specified.");
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
				ThrowSAXException("tauk found before number of G taus specified.");
            tauk=new double[ntaus];
        }
        currentTauk++;
        if(currentTauk>ntaus)
			ThrowSAXException("Too many G taus given.");
        return((char *)&tauk[currentTauk-1]);
    }

    else if(strcmp(xName,"ntausK")==0)
    {   input=INT_NUM;
        return((char *)&ntausK);
    }
    
    else if(strcmp(xName,"Kk")==0)
    {   if(Kk==NULL)
        {   if(ntaus<=0)
                ThrowSAXException("Kk found before number of K taus specified.");
            Kk=new double[ntausK];
        }
        currentKk++;
        if(currentKk>ntausK)
            ThrowSAXException("Too many Kk's given.");
        return UnitsController::ScaledPtr((char *)&Kk[currentKk-1],gScaling,1.e6);
    }
    
    else if(strcmp(xName,"tauKk")==0)
    {   if(tauKk==NULL)
        {   if(ntausK<=0)
                ThrowSAXException("tauKk found before number of K taus specified.");
            tauKk=new double[ntausK];
        }
        currentTauKk++;
        if(currentTauKk>ntausK)
            ThrowSAXException("Too many K taus given.");
        return((char *)&tauKk[currentTauKk-1]);
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
	    
    else if(strcmp(xName,"Tref")==0)
    {   input=DOUBLE_NUM;
        return((char *)&Tref);
    }
        
    else if(strcmp(xName,"C1")==0)
    {   input=DOUBLE_NUM;
        return((char *)&C1base10);
    }
        
    else if(strcmp(xName,"C2")==0)
    {   input=DOUBLE_NUM;
        return((char *)&C2);
    }
    
    else if(strcmp(xName,"mref")==0)
    {   input=DOUBLE_NUM;
        return((char *)&mref);
    }
        
    else if(strcmp(xName,"Cm1")==0)
    {   input=DOUBLE_NUM;
        return((char *)&Cm1base10);
    }
        
    else if(strcmp(xName,"Cm2")==0)
    {   input=DOUBLE_NUM;
        return((char *)&Cm2base10);
    }

    return(MaterialBase::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
// throws std::bad_alloc
const char *Viscoelastic::VerifyAndLoadProperties(int np)
{
    if(pressureLaw!=LINEAR_PRESSURE && pressureLaw!=MGEOS_PRESSURE)
        return "Invalid pressure law was selected";
    
	// Shear modulus - must specify ntaus even if it is zero
    if(currentGk<ntaus || currentTauk<ntaus)
		return "Insufficient Gk or tauk for expected number of G taus.";
    if(ntaus<0)
        return "Number of G taus was never entered.";
    
    // bulk modulus. Not entering ntauKs implies time independent bulk modulus
    if(currentKk<ntausK || currentTauKk<ntausK)
        return "Insufficient Kk or tauKk for expected number of K taus.";
    if(ntausK>0)
    {   if(pressureLaw!=LINEAR_PRESSURE)
            return "Time-dependent bulk modulus requires linear pressure law";
        pressureLaw = TIME_DEPENDENT_PRESSURE;
    }
    else
        ntausK = 0;
    
    cout << "Pressure Law: " << pressureLaw << endl;
    
    // Needs non-zero bulk modulus
    if(K<0)
		return "Required bulk modulus must be positive.";
    
    // zero time shear modulus
    Gered = G0;
	TwoGkred = new double[ntaus];
    for(int k=0;k<ntaus;k++)
    {   Gered += Gk[k];
		TwoGkred[k] = 2.*Gk[k]/rho;
    }
	
	// Convert to specific moduli
	Gered /= rho;

    // bulk modulus
    if(pressureLaw==LINEAR_PRESSURE)
    {    Kered = K/rho;
    }
    
    else if(pressureLaw==TIME_DEPENDENT_PRESSURE)
    {   // zero time shear modulus
        Kered = K;
        Kkred = new double[ntausK];
        for(int k=0;k<ntausK;k++)
        {   Kered += Kk[k];
            Kkred[k] = Kk[k]/rho;
        }
        
        // Convert to specific moduli
        Kered /= rho;
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
			warnExcessiveX = warnings.CreateWarning("Compressive strain has exceeded MGEOS law range in Viscoleastic material",-1,0);
		Xmax = GetMGEOSXmax(gamma0,S1,S2,S3,Kmax);
	}
	
	else
		return "Invalid option for the pressure law";
	
    // to absolute CTE and CME
    CTE = 1.e-6*aI;
    CME = betaI*concSaturation;
    
    // for Cp-Cv (units nJ/(g-K^2))
    Ka2sp = 9.*Kered*CTE*CTE;
    
    // WLF coefficients convert to ln aT
    C1 = log(10.)*C1base10;
    
    // Moisture terms - convert to use ln am = - Cm1*(c-cref)/(Cm2+c) where c = m/mref
    Cm1 = log(10.)*Cm1base10;
    Cm2 = Cm2base10/concSaturation;
    mref /= concSaturation;
    if(mref>=0. && Cm2<=0.)
        return "Cm2 must be greater than zero";

	// call super class
	return MaterialBase::VerifyAndLoadProperties(np);
}

// plane stress not allowed in viscoelasticity
// throws CommonException()
void Viscoelastic::ValidateForUse(int np) const
{	if(np==PLANE_STRESS_MPM)
	{	if(pressureLaw==MGEOS_PRESSURE)
		{	throw CommonException("Viscoelastic materials in plane stress require linear pressure model",
							  		"Viscoelastic::ValidateForUse");
		}
		if(artificialViscosity)
		{	throw CommonException("Viscoelastic materials in plane stress do not support artificial viscosity",
								  "Viscoelastic::ValidateForUse");
		}
	}

	//call super class (why can't call super class?)
	MaterialBase::ValidateForUse(np);
}

#pragma mark Viscoelastic::History Data Methods

// create and return pointer to history variables
// initialize all to zero
// throws std::bad_alloc
char *Viscoelastic::InitHistoryData(char *pchr,MPMBase *mptr)
{
    // count the variables
    numJHistory = pressureLaw==MGEOS_PRESSURE ? 2 : 0;

    if(fmobj->IsThreeD())
        numHistory = numJHistory + 6*ntaus + ntausK;
    else
        numHistory = numJHistory + 4*ntaus + ntausK;
    
    // exit if none
    if(numHistory==0) return NULL;
    
    // all zeros
    double *p = CreateAndZeroDoubles(pchr,numHistory);
    
    // J history starts at 1
    if(numJHistory>0)
    {   p[MGJ_HISTORY] = 1.;
        p[MGJRES_HISTORY] = 1.;
    }
    
    // return pointer
    return (char *)p;
}

// reset history data
void Viscoelastic::ResetHistoryData(char *pchr,MPMBase *mptr)
{	ZeroDoubles(pchr,numHistory);
	// J history starts at 1
	if(numJHistory>0)
	{	double *p = (double *)pchr;
		p[MGJ_HISTORY] = 1.;
		p[MGJRES_HISTORY] = 1.;
	}
}

// Number of history variables - only the plastic law
int Viscoelastic::NumberOfHistoryDoubles(void) const { return numHistory; }

#pragma mark Viscoelastic::Methods

/* Take increments in strain and calculate new Particle: strains, rotation strain,
	stresses, strain energy,
	du are (gradient rates X time increment) to give deformation gradient change
	For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMETRIC_MPM, otherwise dvzz=0
	This material tracks pressure and stores deviatoric stress only in particle stress tensor
 */
void Viscoelastic::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res,int historyOffset) const
{
	// Note: cannot call generic method because need detdF and Vrot below
	
	// previous deformation gradient and stretch
	Matrix3 pFnm1 = mptr->GetDeformationGradientMatrix();
	
	// get incremental deformation gradient and decompose it
	const Matrix3 dF = du.Exponential(incrementalDefGradTerms);
	
	// Update total deformation gradient (saved on particle at the end)
	Matrix3 pF = dF*pFnm1;
	
	// decompose to get previous Rn and Rn-1 and current V
	Matrix3 Rnm1,Rn;
	pFnm1.LeftDecompose(&Rnm1,NULL);
	Matrix3 Vrot = pF.LeftDecompose(&Rn,NULL);
	Matrix3 dR = Rn*Rnm1.Transpose();
	
	// get strain increments in current configuration (dF-dR)F(n-1)Rn^T
	Matrix3 dFmdR = dF - dR;
	Matrix3 detot = dFmdR*(pFnm1*Rn.Transpose());
    
    // shift factors based on mptr->pPreviousTemperature and mptr->pPreviousConcentration
#ifdef TOTAL_STRESS_CALC
    double bshift=1.;
    double pGered = bshift*Gered;
    double pKered = bshift*Kered;
#else
    double bshift=1.;
#endif
	
    // Effective strain by deducting thermal strain (no shear thermal strain because isotropic)
	double eres=CTE*res->dT;
	if(DiffusionTask::HasFluidTransport())
		eres+=CME*res->dC;
	
	// update pressure
	double dispEnergy = 0.,detdF = 1.,J = 1.,dJres = 1.;
	double traceDe = detot.trace();
	double delV = traceDe - 3.*eres;
	
	// history data
	double *ak =(double *)(mptr->GetHistoryPtr(0));
    int ai = numJHistory;

	// find dJ and J if needed (plane stress not allowed)
	if(numJHistory>0)
	{	// large strain volume change
		detdF = dF.determinant();
		J = detdF*ak[MGJ_HISTORY];
		ak[MGJ_HISTORY] = J;
		
		// account for residual strains if needed
		double dJres = exp(3.*eres);
		double Jres = dJres*ak[MGJRES_HISTORY];
		ak[MGJRES_HISTORY] = Jres;
	}
	
	// deviatoric strains increment in de
	// Actually de is finding 2*d(dev e) to avoid many multiplies by two
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
    SubTensor(&ed,&de);

	// increment particle deviatoric stresses - elastic part (factor of 2 in de)
	double dsig[6];
#ifdef TOTAL_STRESS_CALC
	dsig[XX] = pGered*(ed.xx+de.xx);
	dsig[YY] = pGered*(ed.yy+de.yy);
	dsig[ZZ] = pGered*(ed.zz+de.zz);
	dsig[XY] = pGered*(ed.xy+de.xy);
	if(np==THREED_MPM)
	{	dsig[XZ] = pGered*(ed.xz+de.xz);
		dsig[YZ] = pGered*(ed.yz+de.yz);
	}
#else
	dsig[XX] = Gered*de.xx;
	dsig[YY] = Gered*de.yy;
	dsig[ZZ] = Gered*de.zz;
	dsig[XY] = Gered*de.xy;
	if(np==THREED_MPM)
	{	dsig[XZ] = Gered*de.xz;
		dsig[YZ] = Gered*de.yz;
	}
#endif
    
    // get effective time increment
    double delEffTime = GetEffectiveIncrement(mptr,res,delTime,Tref,C1,C2,mref,Cm1,Cm2,0.,0.,0.);
    
	// get internal variable increments, update them, add to incremental stress, and get dissipated energy
    // For plane stress, this gets initial deviatoric stress update only
	Tensor dak;
	int k;
    double omtmp,omtmp2;
    for(k=0;k<ntaus;k++)
    {   GetAlphaArgs(delEffTime,tauk[k],omtmp,omtmp2);
        dak.xx = omtmp*(0.5*ed.xx-ak[ai+XX_HISTORY]) + omtmp2*de.xx;
        dak.yy = omtmp*(0.5*ed.yy-ak[ai+YY_HISTORY]) + omtmp2*de.yy;
        dak.xy = omtmp*(0.5*ed.xy-ak[ai+XY_HISTORY]) + omtmp2*de.xy;
        dak.zz = omtmp*(0.5*ed.zz-ak[ai+ZZ_HISTORY]) + omtmp2*de.zz;
		
		// extra terms for 3D
		if(np==THREED_MPM)
		{	// internal variables
			dak.xz = omtmp*(0.5*ed.xz-ak[ai+XZ_HISTORY]) + omtmp2*de.xz;
			dak.yz = omtmp*(0.5*ed.yz-ak[ai+YZ_HISTORY]) + omtmp2*de.yz;
			
			// update history on particle
			ak[ai+XX_HISTORY] += dak.xx;
			ak[ai+YY_HISTORY] += dak.yy;
            ak[ai+XY_HISTORY] += dak.xy;
			ak[ai+ZZ_HISTORY] += dak.zz;
			ak[ai+XZ_HISTORY] += dak.xz;
			ak[ai+YZ_HISTORY] += dak.yz;
			
			// dissipation (updated deviatoric strain minus updated alpha, remove factor of 2 from strain)
			dispEnergy += bshift*TwoGkred[k]*(dak.xx*(0.5*(ed.xx+de.xx)-ak[ai+XX_HISTORY])
									   + dak.yy*(0.5*(ed.yy+de.yy)-ak[ai+YY_HISTORY])
									   + dak.zz*(0.5*(ed.zz+de.zz)-ak[ai+ZZ_HISTORY])
									   + dak.xy*(0.5*(ed.xy+de.xy)-ak[ai+XY_HISTORY])
									   + dak.xz*(0.5*(ed.xz+de.xz)-ak[ai+XZ_HISTORY])
									   + dak.yz*(0.5*(ed.yz+de.yz)-ak[ai+YZ_HISTORY]));
            
            // next history
            ai += 6;
		}
        else if(np==PLANE_STRESS_MPM)
        {   // history abd dissipated energy done later
            ai += 4;
        }
		else
		{	// update history on particle
			ak[ai+XX_HISTORY] += dak.xx;
			ak[ai+YY_HISTORY] += dak.yy;
            ak[ai+XY_HISTORY] += dak.xy;
			ak[ai+ZZ_HISTORY] += dak.zz;

			// dissipation  (updated deviatoric strain minus updated alpha, remove factor of 2 from strain)
			dispEnergy += bshift*TwoGkred[k]*(dak.xx*(0.5*(ed.xx+de.xx)-ak[ai+XX_HISTORY])
									   + dak.yy*(0.5*(ed.yy+de.yy)-ak[ai+YY_HISTORY])
									   + dak.zz*(0.5*(ed.zz+de.zz)-ak[ai+ZZ_HISTORY])
									   + dak.xy*(0.5*(ed.xy+de.xy)-ak[ai+XY_HISTORY]));
            
            // next history
            ai += 4;
		}
		
#ifdef TOTAL_STRESS_CALC
		// add to stress using updated alphas
		dsig[XX] -= bshift*TwoGkred[k]*ak[ai+XX_HISTORY];
		dsig[YY] -= bshift*TwoGkred[k]*ak[ai+YY_HISTORY];
		dsig[ZZ] -= bshift*TwoGkred[k]*ak[ai+ZZ_HISTORY];
		dsig[XY] -= bshift*TwoGkred[k]*ak[ai+XY_HISTORY];
		if(np==THREED_MPM)
		{	dsig[XZ] -= bshift*TwoGkred[k]*ak[ai+XZ_HISTORY];
			dsig[YZ] -= bshift*TwoGkred[k]*ak[ai+YZ_HISTORY];
		}
#else
		// add to stress increments
		dsig[XX] -= TwoGkred[k]*dak.xx;
		dsig[YY] -= TwoGkred[k]*dak.yy;
		dsig[ZZ] -= TwoGkred[k]*dak.zz;
		dsig[XY] -= TwoGkred[k]*dak.xy;
		if(np==THREED_MPM)
		{	dsig[XZ] -= TwoGkred[k]*dak.xz;
			dsig[YZ] -= TwoGkred[k]*dak.yz;
		}
#endif
	}
    
    double dP=0.,Vstar=0.;
    if(pressureLaw==TIME_DEPENDENT_PRESSURE)
    {   double dTemp=mptr->pPreviousTemperature-thermal.reference;
        double eresStretch=CTE*dTemp;
        if(DiffusionTask::HasDiffusion())
		{   double dConc = diffusion->GetDeltaConcentration(mptr);
            eresStretch = CME*dConc;
        }
    
        // find initial V* (diagonal has 1+eii)
        Vstar = Vrot(0,0)+Vrot(1,1)+Vrot(2,2)-3.-3.*eresStretch;

#ifdef TOTAL_STRESS_CALC
		// elastic pressure
		dP = -bshift*Kered*(Vstar+delV);
#else
        // elastic pressure increment
        dP = -Kered*delV;
#endif
        
        // viscoelastic history
        for(k=0;k<ntaus;k++)
        {   GetAlphaArgs(delEffTime,tauKk[k],omtmp,omtmp2);
            double dalphaV = omtmp*(Vstar-ak[ai]) + 2.*omtmp2*delV;
            
#ifdef TOTAL_STRESS_CALC
            // add to pressure
            dP += bshift*Kkred[k]*(ak[ai]+dalphaV);
#else
			// add to pressure increments
			dP += Kkred[k]*dalphaV;
#endif
            
            if(np!=PLANE_STRESS_MPM)
            {   // update history on particle (plane stress done later)
                ak[ai] += dalphaV;

                // dissipation (updated deviatoric strain minus updated alpha)
                dispEnergy += bshift*Kkred[k]*dalphaV*(Vstar+delV-ak[ai]);
            }
                
            // next history
            ai++;
        }
    }

	// For plane stress, find dezz and adjust all terms
	if(np==PLANE_STRESS_MPM)
	{	double phi = Gered;
        for(k=0;k<ntaus;k++) phi -= bshift*TwoGkred[k]*GetPSArg(delEffTime,tauk[k]);

        double phiK = Kered;
        if(pressureLaw==TIME_DEPENDENT_PRESSURE)
        {   for(k=0;k<ntausK;k++) phiK -= bshift*Kkred[k]*2.*GetPSArg(delEffTime,tauKk[k]);
        }
        else
		{
#ifdef TOTAL_STRESS_CALC
			dP = -bshift*Kered*(Vstar+delV);
#else
			dP = -Kered*delV;
#endif
		}
        
        // dezz
        double dezz = (dP - dsig[ZZ])/(phiK + 4.*phi/3.);
		double thirddezz = dezz/3.;
		
		// adjust deviatoric stress update
		double ds = 2.*phi*thirddezz;
		dsig[XX] -= ds;
		dsig[YY] -= ds;
		dsig[ZZ] += 2.*ds;
		
		// adjust delV for use in pressure update (done below)
		delV += dezz;
		
		// set input strain increment to calculated result (used in work below)
		detot(2,2) = dezz;
		
		// adjust particle deformation gradient (stored below)
		pF(2,2) *= (1.+dezz);
		
		// adjust de for history update and disspated energy (factor of 2 is added)
		de.xx -= 2.*thirddezz;
		de.yy -= 2.*thirddezz;
		de.zz += 4.*thirddezz;
		
		// update history and get dissipation
        ai = numJHistory;
		for(k=0;k<ntaus;k++)
        {   GetAlphaArgs(delEffTime,tauk[k],omtmp,omtmp2);
            dak.xx = omtmp*(0.5*ed.xx-ak[ai+XX_HISTORY]) + omtmp2*de.xx;
            dak.yy = omtmp*(0.5*ed.yy-ak[ai+YY_HISTORY]) + omtmp2*de.yy;
            dak.xy = omtmp*(0.5*ed.xy-ak[ai+XY_HISTORY]) + omtmp2*de.xy;
            dak.zz = omtmp*(0.5*ed.zz-ak[ai+ZZ_HISTORY]) + omtmp2*de.zz;

			// update history on particle
			ak[ai+XX_HISTORY] += dak.xx;
			ak[ai+YY_HISTORY] += dak.yy;
            ak[ai+XY_HISTORY] += dak.xy;
			ak[ai+ZZ_HISTORY] += dak.zz;
			
			// dissipation
			dispEnergy += bshift*TwoGkred[k]*(dak.xx*(0.5*(ed.xx+de.xx)-ak[ai+XX_HISTORY])
									   + dak.yy*(0.5*(ed.yy+de.yy)-ak[ai+YY_HISTORY])
									   + dak.zz*(0.5*(ed.zz+de.zz)-ak[ai+ZZ_HISTORY])
									   + dak.xy*(0.5*(ed.xy+de.xy)-ak[ai+XY_HISTORY]));
            
            // next history
            ai += 4;
		}
        
        if(pressureLaw==TIME_DEPENDENT_PRESSURE)
        {   // change elastic part
            dP -= phiK*dezz;
            
            // update history and get dissipation
            for(k=0;k<ntaus;k++)
            {   GetAlphaArgs(delEffTime,tauKk[k],omtmp,omtmp2);
                double dalphaV = omtmp*(Vstar-ak[ai]) + 2.*omtmp2*delV;
                
                // add to pressure increments
                ak[ai] += dalphaV;
                
                // dissipation (updated volumetric strain minus updated alpha)
                dispEnergy += bshift*Kkred[k]*dalphaV*(Vstar+delV-ak[ai]);
                    
                // next history
                ai++;
            }
        }
	}
	
	// Update particle deviatoric stresses
	Tensor *sp=mptr->GetStressTensor();
	
#ifdef TOTAL_STRESS_CALC
    if(np==THREED_MPM)
    {   // total stress
        Matrix3 stn(dsig[XX],dsig[XY],dsig[XZ],dsig[XY],dsig[YY],dsig[YZ],dsig[XZ],dsig[YZ],dsig[ZZ]);
        
    }
    else
    {   // total stress
        Matrix3 stn(dsig[XX],dsig[XY],dsig[XY],dsig[YY],dsig[ZZ]);
        
    }
#else
	//Tensor st0 = *sp;
	if(np==THREED_MPM)
	{   // incremental rotate of prior stress
		Matrix3 stn(sp->xx,sp->xy,sp->xz,sp->xy,sp->yy,sp->yz,sp->xz,sp->yz,sp->zz);
		Matrix3 str = stn.RMRT(dR);

		if(numJHistory>0)
		{	// convert sigma(n)/rho(n) to sigma(n)/rho(n+1) and add dsigma/rho(n+1)
			sp->xx = detdF*str(0,0)+J*dsig[XX];
			sp->yy = detdF*str(1,1)+J*dsig[YY];
			sp->xy = detdF*str(0,1)+J*dsig[XY];
			sp->zz = detdF*str(2,2)+J*dsig[ZZ];
			sp->yz = detdF*str(1,2)+J*dsig[YZ];
			sp->xz = detdF*str(0,2)+J*dsig[XZ];
		}
		else
		{	// small strain stress increment
			sp->xx = str(0,0)+dsig[XX];
			sp->yy = str(1,1)+dsig[YY];
			sp->xy = str(0,1)+dsig[XY];
			sp->zz = str(2,2)+dsig[ZZ];
			sp->yz = str(1,2)+dsig[YZ];
			sp->xz = str(0,2)+dsig[XZ];
		}
	}
	else
	{	// incremental rotate of prior stress
		Matrix3 stn(sp->xx,sp->xy,sp->xy,sp->yy,sp->zz);
		Matrix3 str = stn.RMRT(dR);
		
		if(numJHistory>0)
		{	// convert sigma(n)/rho(n) to sigma(n)/rho(n+1) and add dsigma/rho(n+1)
			sp->xx = detdF*str(0,0)+J*dsig[XX];
			sp->yy = detdF*str(1,1)+J*dsig[YY];
			sp->xy = detdF*str(0,1)+J*dsig[XY];
			sp->zz = detdF*sp->zz+J*dsig[ZZ];
		}
		else
		{	// small strain stress increment
			sp->xx = str(0,0)+dsig[XX];
			sp->yy = str(1,1)+dsig[YY];
			sp->xy = str(0,1)+dsig[XY];
			sp->zz += dsig[ZZ];
		}
	}
#endif
	
	// incremental work energy = shear energy (dilation and residual energy done in update pressure)
    double shearEnergy = sp->xx*detot(0,0) + sp->yy*detot(1,1) + sp->zz*detot(2,2) + sp->xy*de.xy;
    if(np==THREED_MPM)
    {   shearEnergy += sp->xz*de.xz + sp->yz*de.yz;
    }
    mptr->AddWorkEnergyAndResidualEnergy(shearEnergy,0.);
	
	// finish particle updates
	mptr->SetDeformationGradientMatrix(pF);
	
	// Now update pressure
    double dTq0 = dP;               // pass pressure increment found above, it returns as dTq0
	UpdatePressure(mptr,delV,res,eres,detdF,dJres,delTime,dTq0,dispEnergy);
    
    // dissipated energy per unit mass (dPhi/(rho0 V0)) (Legacy nJ/g)
    mptr->AddPlastEnergy(dispEnergy);
    
    // heat energy
    IncrementHeatEnergy(mptr,dTq0,dispEnergy);
}

// Get parameters to update alpha variables
void Viscoelastic::GetAlphaArgs(double delTime,double tau,double &omtmp,double &omtmp2) const
{   double x = delTime/tau;
    if(x>0.003)
    {   double tmp = exp(-0.5*x);
        omtmp = 1. - tmp*tmp;               // 1-exp(-dt/tau)
        omtmp2 = 0.5*(1. - tmp);            // 0.5*(1-exp(-dt/(2*tau))) because de has factor of 2
    }
    else
    {   // double precision accurate and maybe more precision than 1-exp(-x)
        omtmp = x*(1.-0.5*x*(1-(x/3.)*(1-0.25*x)));         // 1-exp(-dt/tau)
        
        // 0.5*(1-exp(-dt/(2*tau))) because de has factor of 2
        x *= 0.5;
        omtmp2 = 0.5*x*(1.-0.5*x*(1-(x/3.)*(1-0.25*x)));
   }
}

// Get parameter for plane stress calculations
// Get 0.5*(1-exp(-dt/(2*tau))) (0.5 because of 2Gk has factor of 2)
double Viscoelastic::GetPSArg(double delTime,double tau) const
{   double x = 0.5*delTime/tau;
    if(x>0.003)
        return 0.5*(1. - exp(-x));
    
    // double precision accurate and maybe more precision than 1-exp(-x)
    return 0.5*x*(1.-0.5*x*(1.-(x/3.)*(1.-0.25*x)));
}

// This method handles the pressure equation of state. Its tasks are
// 1. Calculate the new pressure
// 2. Update particle pressure
// 3. Increment the particle energy
void Viscoelastic::UpdatePressure(MPMBase *mptr,double delV,ResidualStrains *res,
								  double eres,double detdF,double dJres,
								  double delTime,double &dTq0,double &AVEnergy) const
{
	if(pressureLaw!=MGEOS_PRESSURE)
	{	// pressure change
        double dP = pressureLaw==LINEAR_PRESSURE ? -Kered*delV : dTq0 ;
		
		// artifical viscosity
		// delV is total incremental volumetric strain = total Delta(V)/V
        // not allowed in plane stress because added pressure would change szz
        // may not be valid for time dependent bulk modulus, but is allowed
		if(delV<0. && artificialViscosity)
		{	// Wants K/rho
			double QAVred = GetArtificialViscosity(delV/delTime,sqrt(Kered),mptr);
			AVEnergy += fabs(QAVred*delV);
			dP += QAVred;
		}
		
		// increment pressure
		mptr->IncrementPressure(dP);
		
		// work energy is dU = -P dV + s.de(total)
		// Here do hydrostatic term
		// Work energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
		double avgP = mptr->GetPressure()-0.5*dP;
		mptr->AddWorkEnergyAndResidualEnergy(-avgP*delV,-3.*avgP*eres);
		
		// Isoentropic temperature rise = -(K 3 alpha T)/(rho Cv) (dV/V) (only spot for this material)
		dTq0 = -3.*Kered*CTE*mptr->pPreviousTemperature*(delV+3*eres)/GetHeatCapacity(mptr);
	}
	else
	{	// J is total volume change - may need to reference to free-swelling volume if that works
		// Note that swelling looks like a problem because the sums of strains needs to be adjusted
		//		to stress-free state
		
		// history pointer
		double *h =(double *)(mptr->GetHistoryPtr(0));
		double J = h[MGJ_HISTORY];
		double Jres = h[MGJRES_HISTORY];
		
		// previous pressure
		double P,P0 = mptr->GetPressure();
		double delVMG = 1. - 1./detdF;			// (Vnp1-Vn)/Vnp1
		double Kred;
		
		// M-G EOS
		// Want specific pressure or pressure over current density (using J = rho0/rho)
		if(J<1.)
		{	// new compression J(k+1) = 1-x(k+1)
			double x = 1.-J;
			
			if(x>Xmax)
			{	// law not valid if gets too high
				if(warnings.Issue(warnExcessiveX,-1)==GAVE_WARNING)
				{
#pragma omp critical (output)
					{	cout << "# Excessive x = " << x << " causing Kred to increase more than " << Kmax << " fold" << endl;
						mptr -> Describe();
					}
				}
				
				// reset
				x = Xmax;
			}
			
			// get reduced bulk modulus
			Kred = C0squared*GetMGEOSKRatio(x,gamma0,S1,S2,S3);
			
			// compression law
			// denominator = 1 - S1*x - S2*x^2 - S3*x^3
			double denom = 1./(1. - x*(S1 + x*(S2 + x*S3)));
			
			// current effective and reduced (by rho0) bulk modulus
			double Klawred = C0squared*(1.-0.5*gamma0*x)*denom*denom;
			
			// Pressure from bulk modulus and an energy term
			double e = mptr->GetInternalEnergy();
			P = J*(Klawred*x + gamma0*e);

			// particle isoentropic temperature increment
			dTq0 += -J*gamma0*mptr->pPreviousTemperature*delV;
		}
		else
		{   // In tension hyperelastic law P = - K0(J-1)
			double Jeff = J/Jres;
			Kred = C0squared*Jeff;
			P = -J*C0squared*(Jeff-1.);
			
			// particle isentropic temperature increment
			double Kratio = Jeff;
			dTq0 += -J*Kratio*gamma0*mptr->pPreviousTemperature*delVMG;
		}
		
		// artifical viscosity
		// delVMG is total incremental volumetric strain = total Delta(V)/V
		if(delVMG<0. && artificialViscosity)
		{	double QAVred = GetArtificialViscosity(delVMG/delTime,sqrt(Kred*J),mptr);
			AVEnergy += fabs(QAVred*delVMG);
			P += QAVred;
		}
		
		// set final pressure
		mptr->SetPressure(P);
		
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

// Get effective time increment for viscoelatic materials
double Viscoelastic::GetEffectiveIncrement(MPMBase *mptr,ResidualStrains *res,double dRealTime,
                                           double T0,double cT1,double cT2,double m0,double cC1,
                                           double cC2,double kmsEffect,double cMu1,double cMu2)
{
    // if reference properties not set, using actual time
    if(T0<0. && m0<0.) return dRealTime;
    
    // we want to find R = aT(T)am(T)/(aT(T+dT)am(c+dc))
    // we start ln R = (ln aT(Told)-ln aT(Tnew)) + (ln am(cold)-ln am(cnew)) = lnR
    double lnR = 0.;

    // First check temperature
    bool shifted = false;       // means lnR != 0
    double mlogaold=0.;
    if(T0>=0.)
    {   // previous temperature
        double Tnew = mptr->pPreviousTemperature;
        double Told = Tnew-res->dT;
        
        // freeze viscoleaticity below WLF limit at cT2 degrees below reference temperature
        // May be problem when MS is active
        if(Tnew<=T0-cT2 || Told<=T0-cT2) return 0.;
        
        // only deviates from ln 1=0 only when dT changes
        if(!DbleEqual(res->dT,0.))
        {   lnR = cT1*cT2*res->dT/((cT2+Tnew-T0)*(cT2+Told-T0));
            shifted = true;
        }
        
        // will need this even if does not change (previous -ln aT)
        mlogaold = cT1*(Told-T0)/(cT2+Told-T0);
    }
    
    // Up to here
    // lnR = ln aT(T)/(aT(T+dT)) and mlogaold = - ln aT(T) and shifted=true if lnR!=0

    // now check if moisture effect too (only if solve has a concentration in pDiff[0])
    if(m0>=0. && diffusion!=NULL)
    {   // concentrations from the grid
        double cnew = mptr->pDiff[0]->prevConc;
        double cold = cnew - res->dC;
        
        // only deviates from ln 1=0 and only change ln R when dC changes)
        if(!DbleEqual(res->dC,0.))
        {   // add to temperature changes using (ln am(cold)-ln am(cnew))
            // when done, lnR = ln aT(T)am(c)/(aT(T+dT)am(c+dc))
            double del = cC1*(cC2+m0)*res->dC/((cC2+cnew)*(cC2+cold));
            lnR += del;
            
            // thermal and moisture term
            double scale = 1.;
            double R = exp(lnR);
            if(fabs(R-1.)<0.05)
                scale = (9.+R*(19.-R*(5.-R)))/24.;
            else
                scale = (R-1.)/lnR;
            
            // moisture shift (previous -ln am)
            double mlogamold = cC1*(cold-m0)/(cC2+cold);
            
            // time increment dteff = dt*scale/(aT am)
            double dteff = dRealTime*scale*exp(mlogaold+mlogamold);
        }
        else
        {   // will need this even when dC=0 (previous -ln am); it adds to thermal shift
            mlogaold += cC1*(cold-m0)/(cC2+cold);
        }
    }
    
    // If get here
    // lnR = ln aT(T)am(c)/(aT(T+dT)am(c+dc)) and mlogaold = -ln aT(T)am(c) and shifted=true if lnR!=0
    
    // get increment in effective time (temperature, or temperature and moisture bytut dC=0)
    double scale = 1.;
    if(shifted)
    {   double R = exp(lnR);
        if(fabs(R-1.)<0.05)
            scale = (9.+R*(19.-R*(5.-R)))/24.;
        else
            scale = (R-1.)/lnR;
    }

    // return the increment = dt*scale/(aT am)
    return dRealTime*scale*exp(mlogaold);
}

#pragma mark Viscoelastic::Accessors

// return material type
const char *Viscoelastic::MaterialType(void) const { return "Viscoelastic"; }

// Calculate wave speed in mm/sec
// Uses sqrt((K +4Ge/3)/rho) which is probably the maximum wave speed possible
double Viscoelastic::WaveSpeed(bool threeD,MPMBase *mptr) const { return sqrt((Kered + 4.*Gered/3.)); }

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatory stress into full stress
Tensor Viscoelastic::GetStress(Tensor *sp,double pressure,MPMBase *mptr) const
{	return GetStressPandDev(sp,pressure,mptr);
}

// store a new total stress on a particle's stress and pressure variables
void Viscoelastic::SetStress(Tensor *spnew,MPMBase *mptr) const
{	SetStressPandDev(spnew,mptr);
}

// Increment thickness (zz) stress through deviatoric stress and pressure
void Viscoelastic::IncrementThicknessStress(double dszz,MPMBase *mptr) const
{	IncrementThicknessStressPandDev(dszz,mptr);
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

// Get current relative volume change = J (which this material tracks)
double Viscoelastic::GetCurrentRelativeVolume(MPMBase *mptr,int offset) const
{	if(numJHistory==0) return 1.;
	double *h =(double *)mptr->GetHistoryPtr(offset);
	return h[MGJ_HISTORY];
}

// not supported yet, need to deal with aniostropi properties
bool Viscoelastic::SupportsDiffusion(void) const
{   return DiffusionTask::HasPoroelasticity() ? false : true;
}
