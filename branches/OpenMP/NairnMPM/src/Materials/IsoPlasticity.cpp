/********************************************************************************
    IsoPlasticity.cpp
    NairnMPM
    
    Created by John Nairn, June 16, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		Orthotropic.hpp (TranIsotropic.hpp, MaterialBase.hpp)
********************************************************************************/

#include "IsoPlasticity.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "Materials/LinearHardening.hpp"
#include "Materials/NonlinearHardening.hpp"
#include "Materials/Nonlinear2Hardening.hpp"
#include "Materials/JohnsonCook.hpp"
#include "Materials/SCGLHardening.hpp"
#include "Materials/SLMaterial.hpp"

#pragma mark IsoPlasticity::Constructors and Destructors

// Constructors
IsoPlasticity::IsoPlasticity() {}

// Constructors
IsoPlasticity::IsoPlasticity(char *matName) : IsotropicMat(matName)
{
    plasticLaw = new LinearHardening(this);
}

#pragma mark IsoPlasticity::Initialization

// Read material properties
char *IsoPlasticity::InputMat(char *xName,int &input)
{
    // look for different plastic law
    if(strcmp(xName,"Hardening")==0)
    {	input = HARDENING_LAW_SELECTION;
        return (char *)this;
    }
   
    // check plastic law
    char *ptr = plasticLaw->InputMat(xName,input);
    if(ptr != NULL) return ptr;
    
    // otherwise get material properties
    return(IsotropicMat::InputMat(xName,input));
}

// verify settings and some initial calculations
const char *IsoPlasticity::VerifyAndLoadProperties(int np)
{
	// call plastic law that is used
    const char *ptr = plasticLaw->VerifyAndLoadProperties(np);
    if(ptr != NULL) return ptr;
	
	// check in superclass (along with its initialization)
	ptr = IsotropicMat::VerifyAndLoadProperties(np);
	
	// reduced prooperties
	G0red = C66/rho;
	pr.Gred = G0red;
	pr.Kred = C33/rho - 4.*G0red/3.;							// from C33 = lambda + 2G = K + 4G/3
	
	return ptr;
}

// change hardening law
void IsoPlasticity::SetHardeningLaw(char *lawName)
{   
    // delete old one
    delete plasticLaw;
    plasticLaw = NULL;
    
    // check options
    if(strcmp(lawName,"Linear")==0 || strcmp(lawName,"1")==0)
        plasticLaw = new LinearHardening(this);
    
    else if(strcmp(lawName,"Nonlinear")==0 || strcmp(lawName,"2")==0)
        plasticLaw = new NonlinearHardening(this);
    
    else if(strcmp(lawName,"Nonlinear2")==0 || strcmp(lawName,"6")==0)
        plasticLaw = new Nonlinear2Hardening(this);
    
    else if(strcmp(lawName,"JohnsonCook")==0 || strcmp(lawName,"3")==0)
        plasticLaw = new JohnsonCook(this);
    
    else if(strcmp(lawName,"SCGL")==0 || strcmp(lawName,"4")==0)
        plasticLaw = new SCGLHardening(this);
    
    else if(strcmp(lawName,"SL")==0 || strcmp(lawName,"5")==0)
        plasticLaw = new SLMaterial(this);
    
    // did it work
    if(plasticLaw == NULL)
        ThrowSAXException("The hardening law '%s' is not valid",lawName);
}

// print mechanical properties to the results
void IsoPlasticity::PrintMechanicalProperties(void) const
{	
    IsotropicMat::PrintMechanicalProperties();
	plasticLaw->PrintYieldProperties();
}

// The IsoPlasticity has not history data, but its plasticity law might
char *IsoPlasticity::InitHistoryData(void)
{	int num = plasticLaw->HistoryDoublesNeeded();
	if(num==0) return NULL;
	double *p = CreateAndZeroDoubles(num);
	return (char *)p;
}

#pragma mark IsoPlasticity:Methods

// buffer size for mechanical properties
int IsoPlasticity::SizeOfMechanicalProperties(int &altBufferSize) const
{   altBufferSize = plasticLaw->SizeOfHardeningProps();
    return sizeof(PlasticProperties);
}

// Isotropic material can use read-only initial properties
void *IsoPlasticity::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer) const
{
	PlasticProperties *p = (PlasticProperties *)matBuffer;
	*p = pr;
 	p->hardProps = plasticLaw->GetCopyOfHardeningProps(mptr,np,altBuffer);
	double Gratio = plasticLaw->GetShearRatio(mptr,mptr->GetPressure(),1.,p->hardProps);
	p->Gred = G0red*Gratio;
	
	if(np==PLANE_STRESS_MPM)
	{	// these are terms for plane stress calculations only
		p->psRed = 1./(p->Kred/(2.*p->Gred) + 2./3.);					// (1-2nu)/(1-nu) for plane stress
		p->psLr2G = (p->Kred/(2.*p->Gred) - 1./3.)*p->psRed;			// nu/(1-nu) to find ezz
		p->psKred = p->Kred*p->psRed;									// E/(3(1-v)) to find lambda
	}

	return p;
}

/* Take increments in strain and calculate new Particle: strains, rotation strain, plastic strain,
        stresses, strain energy, plastic energy, dissipated energy, angle
    dvij are (gradient rates X time increment) to give deformation gradient change
   For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMEtRIC_MPM, otherwise dvzz=0
   This is general analysis for isotropic plastic material. Subclass must define
		GetYield() and GetKPrime() and optionally can override more. Those methods
		require history dependent variables and rates (e.g. cum. plastic strain (alpint) and
		plastic strain rate (dalpha/delTime) to be set before they are called.
   This material tracks pressure and stores deviatoric stress only in particle stress
        tensor
*/
void IsoPlasticity::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
    // Effective strain by deducting thermal strain (no shear thermal strain because isotropic)
	//  (note: using unreduced terms in CTE3 and CME3)
	double eres=CTE3*res->dT;
	if(DiffusionTask::active)
		eres+=CME3*res->dC;
	
	// Trial update assuming elastic response
	double delV;
    
    // 3D or 2D
	PlasticProperties *p = (PlasticProperties *)properties;
    if(np==THREED_MPM)
    {   delV = du.trace() - 3.*eres;
        PlasticityConstLaw(mptr,du(0,0),du(1,1),du(2,2),du(0,1),du(1,0),du(0,2),du(2,0),
                           du(1,2),du(2,1),delTime,np,delV,1.,eres,p,res);
        return;
    }
	else if(np==PLANE_STRESS_MPM)
		delV = p->psRed*(du(0,0)+du(1,1)-2.*eres);
	else
		delV = du.trace() - 3.*eres;
    PlasticityConstLaw(mptr,du(0,0),du(1,1),du(0,1),du(1,0),du(2,2),delTime,np,delV,1.,eres,p,res);
}

// To allow some subclasses to support large deformations, the initial calculation for incremental
// deformation gradient (the dvij), volume change (delV) and relative volume (Jnew) can be
// handled first by the subclass. This mtehod then finishes the constitutive law
void IsoPlasticity::PlasticityConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
                                double dvzz,double delTime,int np,double delV,double Jnew,double eres,
									   PlasticProperties *p,ResidualStrains *res) const
{
	// here dvij is total strain increment, dexxr is relative strain by subtracting off eres
    double dexxr = dvxx-eres;			// trial dexx=dvxx
    double deyyr = dvyy-eres;			// trial deyy=dvyy
	double dezzr = dvzz-eres;			// In plane strain trial dezz=0, dvzz=0 except in axisymmetric
    double dgxy = dvxy+dvyx;			// no need to substract residual strain
    double dwrotxy = dvyx-dvxy;

	// allow arbitrary equation of state for pressure
    double P0 = mptr->GetPressure();
	UpdatePressure(mptr,delV,Jnew,np,p,res);
    double Pfinal = mptr->GetPressure();
	
    // Deviatoric stress increment
	Tensor *ep=mptr->GetStrainTensor();
	Tensor *sp=mptr->GetStressTensor();
    Tensor dels,stk,st0=*sp;
    double thirdDelV = delV/3.;
	dels.xx = 2.*p->Gred*(dexxr-thirdDelV);
	dels.yy = 2.*p->Gred*(deyyr-thirdDelV);
	if(np==PLANE_STRESS_MPM)
		dels.zz = Pfinal-P0;
	else
		dels.zz = 2.*p->Gred*(dezzr-thirdDelV);
	dels.xy = p->Gred*dgxy;
    
    // trial deviatoric stress
    stk.xx = st0.xx + dels.xx;
    stk.yy = st0.yy + dels.yy;
	stk.zz = st0.zz + dels.zz;
    stk.xy = st0.xy + dels.xy;
  
    // Calculate plastic potential f = ||s|| - sqrt(2/3)*sy(alpha,rate,...)
	HardeningAlpha alpha;
	plasticLaw->UpdateTrialAlpha(mptr,np,&alpha);			// initialize to last value and zero plastic strain rate
	double strial = GetMagnitudeSFromDev(&stk,np);
    double ftrial = strial - SQRT_TWOTHIRDS*plasticLaw->GetYield(mptr,np,delTime,&alpha,p->hardProps);
	if(ftrial<0.)
	{	// elastic, update stress and strain energy as usual
	
		// Add input strain increment to elastic strain on particle
		ep->xx += dvxx;
		ep->yy += dvyy;
		ep->xy += dgxy;
				
		// increment deviatoric stress (Units N/m^2  cm^3/g) (pressure not needed here)
		Hypo2DCalculations(mptr,-dwrotxy,dels.xx,dels.yy,dels.xy);
		
		// strain energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (uJ/g)
        // energy units are also Pa cm^3/g, i.e., same as stress units
        sp->zz = stk.zz;
		if(np==PLANE_STRAIN_MPM)
        {   mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dexxr + (st0.yy+sp->yy)*deyyr
									   + (st0.xy+sp->xy)*dgxy + (st0.zz+sp->zz)*dezzr));
		}
		else if(np==PLANE_STRESS_MPM)
		{	ep->zz += eres - p->psLr2G*(dexxr+deyyr);
			mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dexxr + (st0.yy+sp->yy)*deyyr
								+ (st0.xy+sp->xy)*dgxy));
		}
		else
		{	ep->zz += dvzz;
			mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dexxr + (st0.yy+sp->yy)*deyyr
									   + (st0.xy+sp->xy)*dgxy + (st0.zz+sp->zz)*dezzr));
		}
		
        // heat energy is Cv(dT-dTq0) - dPhi, but dPhi is zero here
        // and Cv(dT-dTq0) was done in Update pressure
        
		// give subclass material chance to update history variables that change in elastic updates
		ElasticUpdateFinished(mptr,np,delTime);
		
		return;
    }
    
	// Find  lambda for this plastic state
	// Base class finds it numerically, subclass can override if solvable by more efficient methods
    double lambdak = plasticLaw->SolveForLambdaBracketed(mptr,np,strial,&stk,p->Gred,p->psKred,Pfinal,delTime,&alpha,p->hardProps);
	
	// Now have lambda, finish update on this particle
	Tensor dfds;
	if(np==PLANE_STRESS_MPM)
    {   // get final stress state
        double d1 = (1 + p->psKred*lambdak);
		double d2 = (1.+2.*p->Gred*lambdak);
		double n1 = (stk.xx+stk.yy-2.*Pfinal)/d1;
		double n2 = (-stk.xx+stk.yy)/d2;
		double sxx = (n1-n2)/2.;
		double syy = (n1+n2)/2.;
		double txy = stk.xy/d2;
        
        // find increment in deviatoric stress
		dels.xx = sxx+Pfinal-st0.xx;
		dels.yy = syy+Pfinal-st0.yy;
		dels.xy = txy-st0.xy;
        
        // get final direction
        dfds.xx = (2.*sxx-syy)/3.;
        dfds.yy = (2.*syy-sxx)/3.;
        dfds.zz = -(dfds.xx+dfds.yy);
        dfds.xy = txy;				// tensorial shear strain
		
		// If fully plastic, the increment in deviatoric stress should be zero
		// sxx = sigmaxx+Pfinal = st0.xx, syy = sigmayy+Pfinal = st0.yy, szz = Pfinal
		// But first two must be wrong because trace is no longer zero?
		//
		// This is probably wrong
		// i.e.:	n1 + 2*Pfinal = st0.xx+st0.yy
		//			(stk.xx+stk.yy)/d1 - 2*Pfinal*(1/d1-1) = st0.xx+st0.yy
		// But (stk.xx+stk.yy) = (st0.xx+st0.yy) + 2.*Gred*(dexxr+deyyr-2*thirdDelV)
		//			(st0.xx+st0.yy)(1/d1-1) - 2*Pfinal*(1/d1-1) = -2.*Gred*(dexxr+deyyr-2*thirdDelV)/d1
		//			(st0.xx+st0.yy-2*Pfinal)*(1-d1) = -2.*Gred*(dexxr+deyyr-2*thirdDelV)
		//			B*Kred*lam = 2.*Gred*A    or    lam = 2.*Gred*A/(B*Kred)
		// where B = st0.xx+st0.yy-2*Pfinal, and A = (dexxr+deyyr-2*thirdDelV)
		// Then  d1 = 1+2*Gred*A/B, n1 = ((st0.xx+st0.yy) + 2*Gred*A - 2*Pfinal)/d1 = (B+2*Gred*A)*B/(B+2*Gred*A) = B
		//
		// Also expect syy-sxx = sigmayy-sigmaxx = st0.yy-st0.xx = n2
		//			(stk.yy-stk.xx)/d2 = st0.yy-st0.xx
		//			(st0.yy-st0.xx)/d2 + 2*Gred*(deyyr-dexxr)/d2 = (st0.yy-st0.xx)
		//			(st0.yy-st0.xx)(d2-1) = - 2*Gred*(deyyr-dexxr)
		//          lam = (deyyr-dexxr)/(st0.yy-st0.xx)
		// Then n2 = (st0.yy-st0.xx + 2*Gred*(deyyr-dexxr))/(1+2*Gred*(deyyr-dexxr)/(st0.yy-st0.xx)) = st0.yy-st0.xx
		//
		// But these to lam's seem to differ?
		// We can write 
	}
    else
    {   // get final direction
        GetDfDsigma(strial,&stk,np,&dfds);
    }
         
    // Plastic strain increments on particle
    double dexxp = lambdak*dfds.xx;
    double deyyp = lambdak*dfds.yy;
	double dezzp = lambdak*dfds.zz;
    double dgxyp = 2.*lambdak*dfds.xy;     // 2 for engineering plastic shear strain
	Tensor *eplast=mptr->GetPlasticStrainTensor();
    eplast->xx += dexxp;
    eplast->yy += deyyp;
    eplast->xy += dgxyp;
	eplast->zz += dezzp;
    
    // Elastic strain increments on particle
    ep->xx += (dvxx-dexxp);
    ep->yy += (dvyy-deyyp);
    dgxy -= dgxyp;
    ep->xy += dgxy;
	if(np==PLANE_STRESS_MPM)
		ep->zz += eres - p->psLr2G*(dexxr+deyyr+dezzp);
	else
	{	// plain strain and axisymmetric when plastic increments is dezzp
		// In plain strain, dvzz=0 elastic increment is -dezzp to zz strain total is zero
		// Axisymmetry can have non-zero hoop strain
		ep->zz += (dvzz-dezzp);
		dezzr -= dezzp; 
	}
	
    // Elastic strain increment minus the residual terms by now subtracting plastic parts
    dexxr -= dexxp;
    deyyr -= deyyp;
	//dgxy -= dgxyp;			// done above
	//dezzr -= dezzp;			// plane strain and axisymmetry done above, plain stress not needed

	// increment particle deviatoric stresses (plane stress increments found above)
	if(np!=PLANE_STRESS_MPM)
	{	dels.xx -= 2.*p->Gred*dexxp;
		dels.yy -= 2.*p->Gred*deyyp;
		dels.xy -= p->Gred*dgxyp;
		dels.zz -= 2.*p->Gred*dezzp;
		sp->zz += dels.zz;
	}
    else
        sp->zz = Pfinal;          // now equal to Pfinal
        
	Hypo2DCalculations(mptr,-dwrotxy,dels.xx,dels.yy,dels.xy);
	
    // Elastic energy increment per unit mass (dU/(rho0 V0)) (uJ/g)
    double strainEnergy = 0.5*((st0.xx+sp->xx)*dexxr
                               + (st0.yy+sp->yy)*deyyr
                               + (st0.xy+sp->xy)*dgxy);
	if(np!=PLANE_STRESS_MPM)
    {	strainEnergy += 0.5*(st0.zz+sp->zz)*dezzr;
	}
    
    // plastic strain work
    double plastEnergy = lambdak*(sp->xx*dfds.xx + sp->yy*dfds.yy + sp->zz*dfds.zz + 2.*sp->xy*dfds.xy);
    
    // total work
    mptr->AddStrainEnergy(plastEnergy + strainEnergy);
    
    // disispated energy per unit mass (dPhi/(rho0 V0)) (uJ/g)
    double qdalphaTerm = lambdak*SQRT_TWOTHIRDS*plasticLaw->GetYieldIncrement(mptr,np,delTime,&alpha,p->hardProps);
    double dispEnergy = plastEnergy - qdalphaTerm;
    
    // heat energy is Cv(dT-dTq0) - dPhi
    // The dPhi is subtracted here because it will show up in next
    //		time step within Cv dT (if adibatic heating occurs)
    // The Cv(dT-dTq0) was done in update pressure
    IncrementHeatEnergy(mptr,0.,0.,dispEnergy);
    
	// The cumulative dissipated energy is tracked in plastic energy
    // Setting the disp energy allows heating if mechanical energy is on
    mptr->AddPlastEnergy(dispEnergy);
    
	// update internal variables
	plasticLaw->UpdatePlasticInternal(mptr,np,&alpha);
}

// To allow some subclasses to support large deformations, the initial calculation for incremental
// deformation gradient (the dvij), volume change (delV) and relative volume (Jnew) can be
// handled first by the subclass. This mtehod then finishes the constitutive law
void IsoPlasticity::PlasticityConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
        double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np,double delV,double Jnew,double eres,
									   PlasticProperties *p,ResidualStrains *res) const
{
	// here dvij is total strain increment, dexxr is relative strain by subtracting off eres
    double dexxr=dvxx-eres;  
    double deyyr=dvyy-eres;
	double dezzr=dvzz-eres;				// use in plane strain only
    double dgxy=dvxy+dvyx;
    double dgxz=dvxz+dvzx;
	double dgyz=dvyz+dvzy;
	
	// rotational strain increments (particle updated by Hypo3D)
	double dwrotxy=dvyx-dvxy;
	double dwrotxz=dvzx-dvxz;
	double dwrotyz=dvzy-dvyz;
	
	// allow arbitrary equation of state for pressure
	UpdatePressure(mptr,delV,Jnew,np,p,res);
	
    // Elastic deviatoric stress increment
	Tensor *ep=mptr->GetStrainTensor();
	Tensor *sp=mptr->GetStressTensor();
    Tensor stk,st0=*sp;
	double dsig[6];
    double thirdDelV = delV/3.;
	dsig[XX] = 2.*p->Gred*(dexxr-thirdDelV);
	dsig[YY] = 2.*p->Gred*(deyyr-thirdDelV);
	dsig[ZZ] = 2.*p->Gred*(dezzr-thirdDelV);
	dsig[YZ] = p->Gred*dgyz;
	dsig[XZ] = p->Gred*dgxz;
	dsig[XY] = p->Gred*dgxy;
    stk.xx = st0.xx + dsig[XX];
    stk.yy = st0.yy + dsig[YY];
	stk.zz = st0.zz + dsig[ZZ];
    stk.yz = st0.yz + dsig[YZ];
    stk.xz = st0.xz + dsig[XZ];
    stk.xy = st0.xy + dsig[XY];
  
    // Calculate plastic potential f = ||s|| - sqrt(2/3)*sy(alpha,rate,...)
	HardeningAlpha alpha;
	plasticLaw->UpdateTrialAlpha(mptr,np,&alpha);			// initialize to last value
	double strial = GetMagnitudeSFromDev(&stk,np);
    double ftrial = strial - SQRT_TWOTHIRDS*plasticLaw->GetYield(mptr,np,delTime,&alpha,p->hardProps);
	if(ftrial<0.)
	{	// elastic, update stress and strain energy as usual
	
		// Add to total strain
		ep->xx+=dvxx;
		ep->yy+=dvyy;
		ep->zz+=dvzz;
		ep->xy+=dgxy;
		ep->xz+=dgxz;
		ep->yz+=dgyz;
		
		// update stress (need to make hypoelastic)
		Hypo3DCalculations(mptr,dwrotxy,dwrotxz,dwrotyz,dsig);
		
		// strain energy increment per unit mass (dU/(rho0 V0)) (uJ/g)
		mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dexxr
								+ (st0.yy+sp->yy)*deyyr
								+ (st0.zz+sp->zz)*dezzr
								+ (st0.yz+sp->yz)*dgyz
								+ (st0.xz+sp->xz)*dgxz
								+ (st0.xy+sp->xy)*dgxy));
								
        // heat energy is Cv(dT-dTq0) - dPhi, but dPhi is zero here
        // and Cv(dT-dTq0) was done in Update Pressure
        
		// give material chance to update history variables that change in elastic updates
		ElasticUpdateFinished(mptr,np,delTime);
		
		return;
	}
	
	// Find direction of plastic strain and lambda for this plastic state
	// Base class finds it numerically, subclass can override if solvable by more efficient meethods
	Tensor dfds;
	GetDfDsigma(strial,&stk,np,&dfds);
	double lambdak = plasticLaw->SolveForLambdaBracketed(mptr,np,strial,&stk,p->Gred,1.,1.,delTime,&alpha,p->hardProps);
	
	// Now have lambda, finish update on this particle
        
    // Plastic strain increments on particle
    double dexxp = lambdak*dfds.xx;
    double deyyp = lambdak*dfds.yy;
	double dezzp = lambdak*dfds.zz;
    double dgxyp = 2.*lambdak*dfds.xy;     // 2 for engineering plastic shear strain
    double dgxzp = 2.*lambdak*dfds.xz;     // 2 for engineering plastic shear strain
    double dgyzp = 2.*lambdak*dfds.yz;     // 2 for engineering plastic shear strain
	Tensor *eplast=mptr->GetPlasticStrainTensor();
	eplast->xx += dexxp;
    eplast->yy += deyyp;
	eplast->zz += dezzp;
    eplast->xy += dgxyp;
	eplast->xz += dgxzp;
    eplast->yz += dgyzp;
   
    // Elastic strain increments on particle
    ep->xx += (dvxx-dexxp);
    ep->yy += (dvyy-deyyp);
    ep->zz += (dvzz-dezzp);
    dgxy -= dgxyp;
    ep->xy += dgxy;
    dgxz -= dgxzp;
    ep->xz += dgxz;
    dgyz -= dgyzp;
	ep->yz += dgyz;
	
    // Elastic strain increment minus the residual terms by now subtracting plastic parts
    dexxr -= dexxp;
    deyyr -= deyyp;
	dezzr -= dezzp;				// plain strain only
	//dgxy, dgxz, dgyz done above

	// increment particle deviatoric stresses
	dsig[XX] -= 2.*p->Gred*dexxp;
	dsig[YY] -= 2.*p->Gred*deyyp;
	dsig[ZZ] -= 2.*p->Gred*dezzp;
	dsig[YZ] -= p->Gred*dgyzp;
	dsig[XZ] -= p->Gred*dgxzp;
	dsig[XY] -= p->Gred*dgxyp;
	Hypo3DCalculations(mptr,dwrotxy,dwrotxz,dwrotyz,dsig);
	
    // Elastic energy increment per unit mass (dU/(rho0 V0)) (uJ/g)
	double strainEnergy = 0.5*((st0.xx+sp->xx)*dexxr + (st0.yy+sp->yy)*deyyr
                               + (st0.zz+sp->zz)*dezzr + (st0.yz+sp->yz)*dgyz
                               + (st0.xz+sp->xz)*dgxz + (st0.xy+sp->xy)*dgxy);
    
    // plastic strain work
    double plastEnergy = lambdak*(sp->xx*dfds.xx + sp->yy*dfds.yy + sp->zz*dfds.zz
                                  + 2.*sp->xy*dfds.xy + 2.*sp->xz*dfds.xz + 2.*sp->yz*dfds.yz);
    
    // total work
    mptr->AddStrainEnergy(plastEnergy + strainEnergy);
    
    // disispated energy per unit mass (dPhi/(rho0 V0)) (uJ/g)
    double qdalphaTerm = lambdak*SQRT_TWOTHIRDS*plasticLaw->GetYieldIncrement(mptr,np,delTime,&alpha,p->hardProps);
    double dispEnergy = plastEnergy - qdalphaTerm;
    
    // heat energy is Cv(dT-dTq0) - dPhi
    // The dPhi is subtracted here because it will show up in next
    //		time step within Cv dT (if adiabatic heating occurs)
    // The Cv(dT-dTq0) was done already in update pressure
    IncrementHeatEnergy(mptr,0.,0.,dispEnergy);
    
	// The cumulative dissipated energy is tracked in plastic energy
    // Setting the disp energy allows heating if mechanical energy is on
    mptr->AddPlastEnergy(dispEnergy);
    
	// update internal variables
	plasticLaw->UpdatePlasticInternal(mptr,np,&alpha);
}

#pragma mark IsoPlasticity::Custom Methods

// This method handles the pressure equation of state. Its tasks are
// 1. Calculate the new pressure
// 2. Update particle pressure
// 3. Increment the particle energy
// 4. Call plasticLaw to see if it wants to change the shear modulus
// 5. Optionally change delV (which is passed by reference)
// Notes:
//  delV is incremental volume change on this step.
//  J is total volume change at end of step (but it is 1 for low-strain materials)
void IsoPlasticity::UpdatePressure(MPMBase *mptr,double &delV,double J,int np,PlasticProperties *p,ResidualStrains *res) const
{   // pressure change
    double dP = -p->Kred*delV;
    mptr->IncrementPressure(dP);
    
    // work energy is dU = -P dV + s.de(total)
	// Here do hydrostatic term
    // Internal energy increment per unit mass (dU/(rho0 V0)) (uJ/g)
    double avgP = mptr->GetPressure()-0.5*dP;
    mptr->AddStrainEnergy(-avgP*delV);

    // heat energy is Cv dT  - dPhi
	// Here do Cv dT term and dPhi is done later
    IncrementHeatEnergy(mptr,res->dT,0.,0.);
}

// Get magnitude of the deviatoric stress tensor when input is a deviatoric stress
// ||s|| = sqrt(s.s) = sqrt(2J2) where J2 = (1/2)s.s
// In 2D 2J2 = sx^2 + sy^2 + sz^2 + 2*txy^2
// In 3D 2J2 = sx^2 + sy^2 + sz^2 + 2*txy^2 + 2*txz^2 + 2*tyz^2
double IsoPlasticity::GetMagnitudeSFromDev(Tensor *st,int np) const
{
	double s,t;
	
	switch(np)
    {   case THREED_MPM:
            s = st->xx*st->xx + st->yy*st->yy + st->zz*st->zz;
            t = st->xy*st->xy + st->xz*st->xz + st->yz*st->yz;
            break;
            
		default:
			s = st->xx*st->xx + st->yy*st->yy + st->zz*st->zz;
			t = st->xy*st->xy;
			break;
	}
	return sqrt(s+t+t);
}

// return derivatives of the yield function wrt to components of deviatoric stress
// which for isotropic material with f = ||s|| - sqrt(2/3)*sy
// reduces to s/||s|| (written for tensorial plastic strains)
// But do not call for plane stress, which must be found in special case
void IsoPlasticity::GetDfDsigma(double smag,Tensor *st0,int np,Tensor *dfds) const
{
    // s/||s|| = n
    dfds->xx = st0->xx/smag;
    dfds->yy = st0->yy/smag;
    dfds->zz = st0->zz/smag;
    dfds->xy = st0->xy/smag;		// tensorial shear strain
    if(np==THREED_MPM)
    {	dfds->xz = st0->xz/smag;
        dfds->yz = st0->yz/smag;
    }
}

// material can override if history variable changes during elastic update
// (e.g., history dependent on plastic strain rate now zero)
// plastic law might all need it
void IsoPlasticity::ElasticUpdateFinished(MPMBase *mptr,int np,double delTime) const
{	plasticLaw->ElasticUpdateFinished(mptr,np,delTime);
}

#pragma mark IsoPlasticity::Accessors

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatory stress into full stress
Tensor IsoPlasticity::GetStress(Tensor *sp,double pressure) const
{   Tensor stress = *sp;
    stress.xx -= pressure;
    stress.yy -= pressure;
    stress.zz -= pressure;
    return stress;
}

// IsoPlasticity has no history data, by the hardening law might
double IsoPlasticity::GetHistory(int num,char *historyPtr) const
{	return plasticLaw->GetHistory(num,historyPtr);
}

// plastic strain needed to get deformation gradient for this material class
bool IsoPlasticity::PartitionsElasticAndPlasticStrain(void) { return TRUE; }

// Return the material tag
int IsoPlasticity::MaterialTag(void) const { return ISOPLASTICITY; }

// return unique, short name for this material
const char *IsoPlasticity::MaterialType(void) const { return "Isotropic Elastic-Plastic"; }

