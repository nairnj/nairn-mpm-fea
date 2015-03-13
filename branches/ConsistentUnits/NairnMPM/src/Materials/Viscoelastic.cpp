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
}

#pragma mark Viscoelastic::Initialization

// print mechanical properties to output window
void Viscoelastic::PrintMechanicalProperties(void) const
{
    PrintProperty("K",K,"");
	PrintProperty("G0",G0,"");
	PrintProperty("ntaus",(double)ntaus,"");
    cout <<  endl;
	
	int i;
    for(i=0;i<ntaus;i++)
	{	PrintProperty("i",(double)i,"");
		PrintProperty("Gk",Gk[i],"");
		PrintProperty("tauk",1000.*tauk[i],"ms");
        cout << endl;
    }
	
	PrintProperty("a",aI,"");
    cout << endl;
}
    
// Read material properties
char *Viscoelastic::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"G0")==0)
        return((char *)&G0);
    
    else if(strcmp(xName,"K")==0)
        return((char *)&K);
    
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
        return((char *)&Gk[currentGk-1]);
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
		TwoGkred[k] = 2.e6*Gk[k]/rho;
    }
	
	// From MPa to Pa
	Gered *= 1.e6/rho;
	Kered = K*1.e6/rho;
	
	// to absolute CTE and CME
	CTE = 1.e-6*aI;
	CME = betaI*concSaturation;

    // for Cp-Cv (units J/(kg-K^2)
    Ka2sp = 1.e-6*Kered*CTE1*CTE1;
	
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

// create and return pointer to history variables
// initialize all to zero
char *Viscoelastic::InitHistoryData(void)
{
    if(ntaus==0) return NULL;
    
    // allocate array of double pointers (3)
	int blocks;
	if(fmobj->IsThreeD())
		blocks = 6;
	else
		blocks = 4;
    char *p = new char[sizeof(double *)*blocks];
	
    double **h = (double **)p;
    
    // for each allocate ntaus doubles
    //	h[ij_HISTORY][0] to h[ij_HISTORY][ntaus-1] can be read
    //	in MPMConstitutiveLaw() by casting mptr->GetHistoryPtr() pointer as
    //		double **h=(double **)(mptr->GetHistoryPtr())
    h[XX_HISTORY] = new double[ntaus];
    h[YY_HISTORY] = new double[ntaus];
    h[XY_HISTORY] = new double[ntaus];
	h[ZZ_HISTORY] = new double[ntaus];
	if(blocks==6)
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
		if(blocks==6)
		{	h[XZ_HISTORY][k] = 0.;
			h[YZ_HISTORY][k] = 0.;
		}
    }
    
    return p;
}

#pragma mark Viscoelastic::Methods


/* Take increments in strain and calculate new Particle: strains, rotation strain, plastic strain,
		stresses, strain energy, plastic energy, dissipated energy, angle
	dvij are (gradient rates X time increment) to give deformation gradient change
	For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMETRIC_MPM, otherwise dvzz=0
	This material tracks pressure and stores deviatoric stress only in particle stress tensor
*/
void Viscoelastic::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
    // Effective strain by deducting thermal strain (no shear thermal strain because isotropic)
	double eres=CTE*res->dT;
	if(DiffusionTask::active)
		eres+=CME*res->dC;
	
	// update pressure
	double delV = du.trace() - 3.*eres;
	UpdatePressure(mptr,delV,res,eres);
	
	// deviatoric strains increment
	// Actually find 2*de to avoid many multiples by two
	Tensor de;
	double dV = du.trace()/3.;
	de.xx = 2.*(du(0,0) - dV);
	de.yy = 2.*(du(1,1) - dV);
	de.zz = 2.*(du(2,2) - dV);
	de.xy = du(0,1)+du(1,0);
	if(np==THREED_MPM)
	{	de.xz = du(0,2)+du(2,0);
		de.yz = du(1,2)+du(2,1);
	}
	
	// Find initial 2*e(t) in ed
	Tensor *ep=mptr->GetStrainTensor();
	Tensor ed = *ep;
	double thirdV = (ed.xx+ed.yy+ed.zz)/3.;
	ed.xx = 2.*(ed.xx-thirdV);
	ed.yy = 2.*(ed.yy-thirdV);
	ed.zz = 2.*(ed.zz-thirdV);
							 
	// Increment total strains on the particle
    ep->xx += du(0,0);
    ep->yy += du(1,1);
	ep->zz += du(2,2);
    ep->xy += de.xy;
	if(np==THREED_MPM)
	{	ep->xz += de.xz;
		ep->yz += de.yz;
	}
	
	// get internal variable increments and update them too
	Tensor dak;
    double **ak =(double **)(mptr->GetHistoryPtr());
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
		
		if(np==THREED_MPM)
		{	dak.xz = tmpm1*ak[XZ_HISTORY][k] + arg*(tmpp1*ed.xz + de.xz);
			dak.yz = tmpm1*ak[YZ_HISTORY][k] + arg*(tmpp1*ed.yz + de.yz);
			ak[XZ_HISTORY][k] += dak.xz;
			ak[YZ_HISTORY][k] += dak.yz;
		}
    }
	
	// increment particle deviatoric stresses
	double dsig[6];
	dsig[XX] = Gered*de.xx;
	dsig[YY] = Gered*de.yy;
	dsig[ZZ] = Gered*de.zz;
	dsig[XY] = Gered*de.xy;
	if(np==THREED_MPM)
	{	dsig[XZ] = Gered*de.xz;
		dsig[YZ] = Gered*de.yz;
	}
    for(k=0;k<ntaus;k++)
	{	dsig[XX] -= TwoGkred[k]*dak.xx;
		dsig[YY] -= TwoGkred[k]*dak.yy;
		dsig[ZZ] -= TwoGkred[k]*dak.zz;
		dsig[XY] -= TwoGkred[k]*dak.xy;
		if(np==THREED_MPM)
		{	dsig[XZ] -= TwoGkred[k]*dak.xz;
			dsig[YZ] -= TwoGkred[k]*dak.yz;
		}
	}
	
	// Hypoelastic increment of particle deviatoric stresses
	Tensor *sp=mptr->GetStressTensor();
	Tensor st0 = *sp;
	double dwrotxy = du(1,0)-du(0,1);
	if(np==THREED_MPM)
	{	double dwrotxz = du(2,0)-du(0,2);
		double dwrotyz = du(2,1)-du(1,2);
		Hypo3DCalculations(mptr,dwrotxy,dwrotxz,dwrotyz,dsig);
	}
	else
	{	Hypo2DCalculations(mptr,-dwrotxy,dsig[XX],dsig[YY],dsig[XY]);
		sp->zz += dsig[ZZ];
	}
	
	// incremental work energy = shear energy (dilation and residual energy done in update pressure)
    double shearEnergy = 0.5*((sp->xx+st0.xx)*du(0,0) + (sp->yy+st0.yy)*du(1,1) + (sp->zz+st0.zz)*du(2,2)+
							  (sp->xy+st0.xy)*de.xy);
    if(np==THREED_MPM)
    {   shearEnergy += 0.5*((sp->xz+st0.xz)*de.xz + (sp->yz+st0.yz)*de.yz);
    }
    mptr->AddWorkEnergyAndResidualEnergy(shearEnergy,0.);
	
    // disispated energy per unit mass (dPhi/(rho0 V0)) (uJ/g)
    double dispEnergy = 0.;
    for(k=0;k<ntaus;k++)
	{	dispEnergy += TwoGkred[k]*(dak.xx*(0.5*(ed.xx+0.5*de.xx)-ak[XX_HISTORY][k]+0.5*dak.xx)
						+ dak.yy*(0.5*(ed.yy+0.5*de.yy)-ak[YY_HISTORY][k]+0.5*dak.yy)
						+ dak.zz*(0.5*(ed.zz+0.5*de.zz)-ak[ZZ_HISTORY][k]+0.5*dak.zz)
						+ dak.xy*(0.5*(ed.xy+0.5*de.xy)-ak[XY_HISTORY][k]+0.5*dak.xy));
		if(np==THREED_MPM)
		{	dispEnergy += TwoGkred[k]*(dak.xz*(0.5*(ed.xz+0.5*de.xz)-ak[XZ_HISTORY][k]+0.5*dak.xz)
							+ dak.yz*(0.5*(ed.yz+0.5*de.yz)-ak[YZ_HISTORY][k]+0.5*dak.yz));
		}
	}
    mptr->AddPlastEnergy(dispEnergy);
    
    // heat energy is Cv dT - dPhi
    // The dPhi is subtracted here because it will show up in next
    //		time step within Cv dT (if adibatic heating occurs)
    // The Cv dT was done in update pressure
    IncrementHeatEnergy(mptr,0.,0.,dispEnergy);
}

// This method handles the pressure equation of state. Its tasks are
// 1. Calculate the new pressure
// 2. Update particle pressure
// 3. Increment the particle energy
void Viscoelastic::UpdatePressure(MPMBase *mptr,double &delV,ResidualStrains *res,double eres) const
{   // pressure change
    double dP = -Kered*delV;
    mptr->IncrementPressure(dP);
    
    // work energy is dU = -P dV + s.de(total)
	// Here do hydrostatic term
    // Work energy increment per unit mass (dU/(rho0 V0)) (uJ/g)
    double avgP = mptr->GetPressure()-0.5*dP;
    mptr->AddWorkEnergyAndResidualEnergy(-avgP*delV,-3.*avgP*eres);
	
    // heat energy is Cv dT  - dPhi
	// Here do Cv dT term and dPhi is done later
    IncrementHeatEnergy(mptr,res->dT,0.,0.);
}

// convert J to K using isotropic method
Vector Viscoelastic::ConvertJToK(Vector d,Vector C,Vector J0,int np)
{	double nuLS = (3.*Kered-2.*Gered)/(6.*Kered+2.*Gered);
	return IsotropicJToK(d,C,J0,np,nuLS,Gered*rho*1.e-6);
}

// From thermodyanamics Cp-Cv = 9 K a^2 T/rho
// Units mJ/(g-K) = J/(kg-m)
// Here using K0 and rho0 - could modify if needed
double Viscoelastic::GetCpMinusCv(MPMBase *mptr) const
{   return mptr!=NULL ? Ka2sp*mptr->pPreviousTemperature : Ka2sp*thermal.reference;
}

#pragma mark Viscoelastic::Accessors

// return material type
const char *Viscoelastic::MaterialType(void) const { return "Viscoelastic"; }

// Return the material tag
int Viscoelastic::MaterialTag(void) const { return VISCOELASTIC; }

/* Calculate wave speed in mm/sec (because G in MPa and rho in g/cm^3)
	Uses sqrt((K +4Ge/3)/rho) which is probably the maximum wave speed possible
*/
double Viscoelastic::WaveSpeed(bool threeD,MPMBase *mptr) const { return sqrt(1.e3*(Kered + 4.*Gered/3.)); }

// Should support archiving history - if it is useful
double Viscoelastic::GetHistory(int num,char *historyPtr) const { return (double)0; }

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatory stress into full stress
Tensor Viscoelastic::GetStress(Tensor *sp,double pressure) const
{
	Tensor stress = *sp;
    stress.xx -= pressure;
    stress.yy -= pressure;
    stress.zz -= pressure;
    return stress;
}

