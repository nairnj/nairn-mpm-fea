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

#pragma mark IsoPlasticity::Constructors and Destructors

// Constructors
IsoPlasticity::IsoPlasticity() {}

// Constructors
IsoPlasticity::IsoPlasticity(char *matName) : IsotropicMat(matName)
{
	readYield=FALSE;
}

#pragma mark IsoPlasticity::Initialization

// Read material properties
char *IsoPlasticity::InputMat(char *xName,int &input)
{
    if(strcmp(xName,"yield")==0)
    {	input=DOUBLE_NUM;
        readYield=TRUE;
        return((char *)&yield);
    }
    
    return(IsotropicMat::InputMat(xName,input));
}

// verify settings and some initial calculations
const char *IsoPlasticity::VerifyProperties(int np)
{	
	// check properties
	if(!readYield) return "The yield stress is missing";

	// call super class
	return IsotropicMat::VerifyProperties(np);
}

// Private properties used in constitutive law
// For variable shear and bulk moduli, subclass can overrive
//		LoadMechanicalProps(MPMBase *mptr,int np) and set new
//		Gred and Kred
// Here gets yldred, Gred, Kred, psRed, psLr2G, and psKred
void IsoPlasticity::InitialLoadMechProps(int makeSpecific,int np)
{
	// reduced prooperties
    yldred = yield*1.e6/rho;
	Gred = C66/rho;
	Kred = C33/rho - 4.*Gred/3.;	// from C33 = lambda + 2G = K + 4G/3
	
	// these are terms for plane stress calculations only
	psRed = 1./(Kred/(2.*Gred) + 2./3.);			// (1-2nu)/(1-nu) for plane stress
	psLr2G = (Kred/(2.*Gred) - 1./3.)*psRed;		// nu/(1-nu) to find ezz
	psKred = Kred*psRed;							// E/(3(1-v)) to find lambda
	
	// nothing needed from superclasses
}

// The base class history variable is cummulative equivalent plastic strain
//		(defined as dalpha = sqrt((2/3)||dep||))
// If super class needs to override this method, always save the first double
//		for the cumulative plastic strain.
char *IsoPlasticity::MaterialData(void)
{
	double *p=new double;
	*p=0.;
	return (char *)p;
}

#pragma mark IsoPlasticity:Methods

/* For 2D MPM analysis, take increments in strain and calculate new
    Particle: strains, rotation strain, plastic strain, stresses, strain energy, 
		plastic energy, dissipated energy, angle
    dvij are (gradient rates X time increment) to give deformation gradient change
	This is general analysis for isotropic plastic material. Subclass must define
		GetYield() and GetKPrime() and optionally can override more. Those methods
		require history dependent variables and rates (e.g. cum. plastic strain (alpint) and
		plastic strain rate (dalpha/delTime) to be set before they are called.
*/
void IsoPlasticity::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
        double delTime,int np)
{
    // Effective strain by deducting thermal strain (no shear thermal strain because isotropic)
	//  (note: using unreduced terms in CTE3 and CME3)
	double eres=CTE3*ConductionTask::dTemperature;
	if(DiffusionTask::active)
		eres+=CME3*DiffusionTask::dConcentration;
	
	// here dvij is total strain increment, dexxr is relative strain by subtracting off eres
    double dexxr=dvxx-eres;			// trial dexx=dvxx
    double deyyr=dvyy-eres;			// trial deyy=dvyy
	double dezzr=-eres;				// used in plane strain only, where trial dezz=0
    double dgxy=dvxy+dvyx;			// no need to substract residual strain
	
	// Trial update assuming elastic response
	double delV,dP;
	if(np==PLANE_STRAIN_MPM)
		delV = (dexxr+deyyr+dezzr)/3.;			// = (dvxx+dvyy)/3. - eres
	else
		delV = psRed*(dexxr+deyyr)/3.;
	
	// allow arbitrary equation of state for pressure
	dP = GetPressureChange(mptr,delV,np);
	
    // Elastic stress increment
	Tensor *ep=mptr->GetStrainTensor();
	Tensor *sp=mptr->GetStressTensor();
    Tensor dels,stk,st0=*sp;
	dels.xx = 2.*Gred*(dexxr-delV) - dP;
	dels.yy = 2.*Gred*(deyyr-delV) - dP;
	if(np==PLANE_STRAIN_MPM)
		dels.zz = 2.*Gred*(dezzr-delV) - dP;
	else
		dels.zz = 0.;
	dels.xy = Gred*dgxy;
    stk.xx = st0.xx + dels.xx;
    stk.yy = st0.yy + dels.yy;
	stk.zz = st0.zz + dels.zz;
    stk.xy = st0.xy + dels.xy;
  
    // Calculate plastic potential f = ||s|| - sqrt(2/3)*sy(alpha,rate,...)
	UpdateTrialAlpha(mptr,np);			// initialize to last value and zero plastic strain rate
	double strial = GetMagnitudeS(&stk,np);
	double ftrial = strial - SQRT_TWOTHIRDS*GetYield(mptr,np,delTime);
	if(ftrial<0)
	{	// elastic, update stress and strain energy as usual
	
		// Add input strain increment to elastic strain on particle
		ep->xx += dvxx;
		ep->yy += dvyy;
		ep->xy += dgxy;
		double dwrotxy=dvyx-dvxy;
				
		// increment stress (Units N/m^2  cm^3/g)
		Hypo2DCalculations(mptr,-dwrotxy,dels.xx,dels.yy,dels.xy);
		
		// out of plane
		if(np==PLANE_STRAIN_MPM)
			sp->zz=stk.zz;
		else
			ep->zz += eres - psLr2G*(dexxr+deyyr);
		
		// strain energy (by midpoint rule)
		mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dexxr + (st0.yy+sp->yy)*deyyr
								+ (st0.xy+sp->xy)*dgxy));
		if(np==PLANE_STRAIN_MPM)
			mptr->AddStrainEnergy(0.5*(st0.zz+sp->zz)*dezzr);
		
		// give material chance to update history variables that change in elastic updates
		ElasticUpdateFinished(mptr,np,delTime);
		
		return;
    }
	
	// Find  lambda for this plastic state
	// Base class finds it numerically, subclass can override if solvable by more efficient methods
    double lambdak = SolveForLambdaBracketed(mptr,np,strial,&stk,delTime);
	
	// Now have lambda, finish update on this particle
	if(np==PLANE_STRESS_MPM)
	{	double d1 = (1 + psKred*lambdak);
		double d2 = (1.+2.*Gred*lambdak);
		double n1 = (stk.xx+stk.yy)/d1;
		double n2 = (-stk.xx+stk.yy)/d2;
		stk.xx = (n1-n2)/2.;
		stk.yy = (n1+n2)/2.;
		stk.xy = stk.xy/d2;
		UpdateTrialAlpha(mptr,np,lambdak,GetMagnitudeS(&stk,np));
		dels.xx=stk.xx-st0.xx;
		dels.yy=stk.yy-st0.yy;
		dels.xy=stk.xy-st0.xy;
	}
	
	// get final direction
	GetDfDsigma(strial,&stk,np);
         
    // Plastic strain increments on particle
    double dexxp = lambdak*dfdsxx;
    double deyyp = lambdak*dfdsyy;
	double dezzp = lambdak*dfdszz;
    double dgxyp = 2.*lambdak*dfdtxy;     // 2 for engineering plastic shear strain
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
		ep->zz += eres - psLr2G*(dexxr+deyyr+dezzp);
	else
	{	ep->zz -= dezzp;		// += (0-dezzp) where 0 is the trial dezz
		dezzr -= dezzp;
	}
	
	// rotational strain
	double dwrotxy=dvyx-dvxy;
	
    // Elastic strain increment minus the residual terms by now subtracting plastic parts
    dexxr -= dexxp;
    deyyr -= deyyp;
	//dgxy -= dgxyp;			// done above
	//dezzr -= dezzp;			// plane strain only done above

	// increment particle stresses (plane stress increments found above)
	// First since plastic strain does not change volume or pressure, we can
	//    use the trial delV and dP from above
	if(np==PLANE_STRAIN_MPM)
	{	dels.xx -= 2.*Gred*dexxp;
		dels.yy -= 2.*Gred*deyyp;
		dels.xy -= Gred*dgxyp;
		dels.zz -= 2.*Gred*dezzp;
		sp->zz += dels.zz;
	}
	Hypo2DCalculations(mptr,-dwrotxy,dels.xx,dels.yy,dels.xy);
	
    // Elastic energy density
    mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dexxr
                        + (st0.yy+sp->yy)*deyyr
                        + (st0.xy+sp->xy)*dgxy));

    // Plastic energy increment
	double dispEnergy=0.5*((st0.xx+sp->xx)*dexxp
                        + (st0.yy+sp->yy)*deyyp
                        + (st0.xy+sp->xy)*dgxyp);
	
	if(np==PLANE_STRAIN_MPM)
    {	mptr->AddStrainEnergy(0.5*(st0.zz+sp->zz)*dezzr);
		dispEnergy += 0.5*(st0.zz+sp->zz)*dezzp;
	}

	// add plastic energy to the particle
	mptr->AddDispEnergy(dispEnergy);
    mptr->AddPlastEnergy(dispEnergy);
	
	// update internal variables
	UpdatePlasticInternal(mptr,np);
}

/* For 3D MPM analysis, take increments in strain and calculate new
    Particle: strains, rotation strain, stresses, strain energy, angle
    dvij are (gradient rates X time increment) to give deformation gradient change
    Assumes linear elastic, uses hypoelastic correction
*/
void IsoPlasticity::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
        double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np)
{
    // Effective strain by deducting thermal strain (no shear thermal strain because isotropic)
	//  (note: using unreduced terms in CTE3 and CME3)
	double eres=CTE3*ConductionTask::dTemperature;
	if(DiffusionTask::active)
		eres+=CME3*DiffusionTask::dConcentration;
	
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
	
	// Trial update assuming elastic response
	double delV = (dexxr+deyyr+dezzr)/3.;			// = (dvxx+dvyy+dvyy)/3. - eres
	
	// allow arbitrary equation of state for pressure
	double dP = GetPressureChange(mptr,delV,np);
	
    // Elastic stress increment
	Tensor *ep=mptr->GetStrainTensor();
	Tensor *sp=mptr->GetStressTensor();
    Tensor stk,st0=*sp;
	double dsig[6];
	dsig[XX] = 2.*Gred*(dexxr-delV) - dP;
	dsig[YY] = 2.*Gred*(deyyr-delV) - dP;
	dsig[ZZ] = 2.*Gred*(dezzr-delV) - dP;
	dsig[YZ] = Gred*dgyz;
	dsig[XZ] = Gred*dgxz;
	dsig[XY] = Gred*dgxy;
    stk.xx = st0.xx + dsig[XX];
    stk.yy = st0.yy + dsig[YY];
	stk.zz = st0.zz + dsig[ZZ];
    stk.yz = st0.yz + dsig[YZ];
    stk.xz = st0.xz + dsig[XZ];
    stk.xy = st0.xy + dsig[XY];
  
    // Calculate plastic potential f = ||s|| - sqrt(2/3)*sy(alpha,rate,...)
	UpdateTrialAlpha(mptr,np,(double)0.,(double)1.);			// initialize to last value
	double strial = GetMagnitudeS(&stk,np);
	double ftrial = strial - SQRT_TWOTHIRDS*GetYield(mptr,np,delTime);
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
		
		// strain energy
		mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dexxr
								+ (st0.yy+sp->yy)*deyyr
								+ (st0.zz+sp->zz)*dezzr
								+ (st0.yz+sp->yz)*dgyz
								+ (st0.xz+sp->xz)*dgxz
								+ (st0.xy+sp->xy)*dgxy));
								
		// give material chance to update history variables that change in elastic updates
		ElasticUpdateFinished(mptr,np,delTime);
		
		return;
	}
	
	// Find direction of plastic strain and lambda for this plastic state
	// Base class finds it numerically, subclass can override if solvable by more efficient meethods
	GetDfDsigma(strial,&stk,np);
	double lambdak = SolveForLambdaBracketed(mptr,np,strial,&stk,delTime);
	
	// Now have lambda, finish update on this particle
        
    // Plastic strain increments on particle
    double dexxp = lambdak*dfdsxx;
    double deyyp = lambdak*dfdsyy;
	double dezzp = lambdak*dfdszz;
    double dgxyp = 2.*lambdak*dfdtxy;     // 2 for engineering plastic shear strain
    double dgxzp = 2.*lambdak*dfdtxz;     // 2 for engineering plastic shear strain
    double dgyzp = 2.*lambdak*dfdtyz;     // 2 for engineering plastic shear strain
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

	// increment particle stresses
	dsig[XX] -= 2.*Gred*dexxp;
	dsig[YY] -= 2.*Gred*deyyp;
	dsig[ZZ] -= 2.*Gred*dezzp;
	dsig[YZ] -= Gred*dgyzp;
	dsig[XZ] -= Gred*dgxzp;
	dsig[XY] -= Gred*dgxyp;
	Hypo3DCalculations(mptr,dwrotxy,dwrotxz,dwrotyz,dsig);
	
    // Elastic energy density
	mptr->AddStrainEnergy(0.5*((st0.xx+sp->xx)*dexxr
							+ (st0.yy+sp->yy)*deyyr
							+ (st0.zz+sp->zz)*dezzr
							+ (st0.yz+sp->yz)*dgyz
							+ (st0.xz+sp->xz)*dgxz
							+ (st0.xy+sp->xy)*dgxy));

    // Plastic energy increment
	double dispEnergy=0.5*(0.5*((st0.xx+sp->xx)*dexxp
							+ (st0.yy+sp->yy)*deyyp
							+ (st0.zz+sp->zz)*dezzp
							+ (st0.yz+sp->yz)*dgyzp
							+ (st0.xz+sp->xz)*dgxzp
							+ (st0.xy+sp->xy)*dgxyp));

	// add plastic energy to the particle
	mptr->AddDispEnergy(dispEnergy);
    mptr->AddPlastEnergy(dispEnergy);
	
	// update internal variables
	UpdatePlasticInternal(mptr,np);
}

#pragma mark IsoPlasticity::Custom Methods

// implement equation of state if needed for this material
double IsoPlasticity::GetPressureChange(MPMBase *mptr,double &delV,int np)
{	return -3.*Kred*delV;
}

/* Get internal variables while iterating to find lambda
	Here alpha (in alpint) is cumulative equivalent plastic strain and
	The increment in alpha is dalpha. dalpha/delTime is plastic strain rate
	For plane strain and 3D
		dalpha = lambda sqrt(2/3) df/dsigma::df/dsigma = lambda sqrt(2/3) since df/dsigma is unit vector
	For plane stress
		dalpha = lambda sqrt(2/3) fnp1
	First one is called at start to initialize alpint and dalpha
*/
void IsoPlasticity::UpdateTrialAlpha(MPMBase *mptr,int np)
{	alpint = mptr->GetHistoryDble();
	dalpha = 0.;
}
void IsoPlasticity::UpdateTrialAlpha(MPMBase *mptr,int np,double lambdak,double fnp1)
{	dalpha = (np==PLANE_STRESS_MPM) ? SQRT_TWOTHIRDS*lambdak*fnp1 : SQRT_TWOTHIRDS*lambdak ;
	alpint = mptr->GetHistoryDble() + dalpha;
}

// The prior soltution for lambda tracked the internal variable. Just move to the particle now
void IsoPlasticity::UpdatePlasticInternal(MPMBase *mptr,int np)
{	mptr->SetHistoryDble(alpint);
}

// Get magnitude of the deviatoric stress tensor
// ||s|| = sqrt(2J2) where J2 = (1/6)((sx-sy)^2+(sy-sz)^2+(sx-sz)^2) + txy^2+txz^2+tyz^2
// In 2D 2J2 = (2/3)(sx^2 + sy^2 - sx*sy + sz*(sz-sx-sy)) + 2*txy^2
// In 3D 2J2 = (2/3)(sx^2 + sy^2 - sx*sy + sz*(sz-sx-sy)) + 2*txy^2 + 2*txz^2 + 2*tyz^2
double IsoPlasticity::GetMagnitudeS(Tensor *st,int np)
{
	double s,t;
	
	switch(np)
	{	case PLANE_STRAIN_MPM:
			s = fmax(st->xx*st->xx + st->yy*st->yy - st->xx*st->yy + st->zz*(st->zz-st->xx-st->yy),0.);
			t = st->xy*st->xy;
			break;
		case PLANE_STRESS_MPM:
			s = fmax(st->xx*st->xx + st->yy*st->yy - st->xx*st->yy,0.);
			t = st->xy*st->xy;
			break;
		case THREED_MPM:
		default:
			s = fmax(st->xx*st->xx + st->yy*st->yy - st->xx*st->yy + st->zz*(st->zz-st->xx-st->yy),0.);
			t = st->xy*st->xy + st->xz*st->xz + st->yz*st->yz;
			break;
	}
	return sqrt(2.*(s/3. + t));
}

// return derivatives of the yield function wrt to components of stress
// which for isotropic material with f = ||s|| - sqrt(2/3)*sy
// reduces to s/||s|| (written for tensorial plastic strains)
void IsoPlasticity::GetDfDsigma(double smag,Tensor *st0,int np)
{
	double p = (st0->xx + st0->yy + st0->zz)/3.;
	if(np==PLANE_STRESS_MPM)
	{	dfdsxx=(st0->xx - p);
		dfdsyy=(st0->yy - p);
		dfdszz=-(dfdsxx+dfdsyy);
		dfdtxy=st0->xy;				// tensorial shear strain
	}
	else
	{	dfdsxx=(st0->xx - p)/smag;
		dfdsyy=(st0->yy - p)/smag;
		dfdszz=(st0->zz - p)/smag;
		dfdtxy=st0->xy/smag;		// tensorial shear strain
		if(np==THREED_MPM)
		{	dfdtxz=st0->xz/smag;
			dfdtyz=st0->yz/smag;
		}
	}
}

/* Solve numerically for lambda
    This method is not used by the IsoPlasticity class, but it may be used by subclasses
        by overriding SolveForLambdaBracketed() and calling this insteas
    Uses Newton's law with initial guess being lambda = dalpha/sqrt(2/3). The solution is not
        bracket, which means it may not be safe. It faster than bracketing when it is safe, but'
        otherwise should not be used (currently used by SLMaterial and VonMisesHardening)
    Set alpint and dalpha before calling
*/
double IsoPlasticity::SolveForLambda(MPMBase *mptr,int np,double strial,Tensor *stk,double delTime)
{
	// initial lambdk from dalpha set before call, often 0, but might be otherwise
	double lambdak=dalpha/SQRT_TWOTHIRDS;
	int step=1;
	
	if(np==PLANE_STRESS_MPM)
	{	double n2trial = -stk->xx+stk->yy;
		n2trial *= n2trial/2;
		n2trial += 2.*stk->xy*stk->xy;
		double n1trial = stk->xx+stk->yy;
		n1trial *= n1trial/6.;
		while(true)
		{	// update iterative variables (lambda, alpha, stress)
			double d1 = (1 + psKred*lambdak);
			double d2 = (1.+2.*Gred*lambdak);
			double fnp12 = n1trial/(d1*d1) + n2trial/(d2*d2);
			double kyld = GetYield(mptr,np,delTime);
			double glam = 0.5*fnp12 - kyld*kyld/3.;
			double fnp1 = sqrt(fnp12);
			double slope = -(psKred*n1trial/(d1*d1*d1) + 2*Gred*n2trial/(d2*d2*d2)) - GetK2Prime(mptr,fnp1,delTime);
			double delLam = -glam/slope;
			lambdak += delLam;
			UpdateTrialAlpha(mptr,np,lambdak,fnp1);
			
			// check for convergence
			if(LambdaConverged(step++,lambdak,delLam)) break;
		}
	}
	else
	{	while(true)
		{	// update iterative variables (lambda, alpha)
			double glam = -SQRT_TWOTHIRDS*GetYield(mptr,np,delTime) + strial - 2*Gred*lambdak;
			double slope = -2.*Gred - GetKPrime(mptr,np,delTime);
			double delLam = -glam/slope;
			lambdak += delLam;
			UpdateTrialAlpha(mptr,np,lambdak,(double)0.);
 			
			// check for convergence
			if(LambdaConverged(step++,lambdak,delLam)) break;
		}
	}
	return lambdak;
}

/* Solve numerically for lambda by safe Newton's method (i.e., with bracketing)
    Subclasses can override for analytical solution is possible or if more efficient method
        is available (e.g., non-bracketed method in SolveForLambda())
    the input ftrial is f function when lambda=0 (but not useful in in plane stress)
*/
double IsoPlasticity::SolveForLambdaBracketed(MPMBase *mptr,int np,double strial,Tensor *stk,double delTime)
{
    double xl,xh;
    BracketSolution(mptr,np,strial,stk,delTime,&xl,&xh);
        
	// initial lambdk midpoint of the brackets
	double lambdak=0.5*(xl+xh);
    UpdateTrialAlpha(mptr,np,lambdak,(double)0.);
    double dxold=fabs(xh-xl);
    double dx=dxold;
	int step=1;
	
	if(np==PLANE_STRESS_MPM)
	{	double n2trial = -stk->xx+stk->yy;
		n2trial *= n2trial/2;
		n2trial += 2.*stk->xy*stk->xy;
		double n1trial = stk->xx+stk->yy;
		n1trial *= n1trial/6.;
        while(true)
        {	// update iterative variables (lambda, alpha)
			double d1 = (1 + psKred*lambdak);
			double d2 = (1.+2.*Gred*lambdak);
			double fnp12 = n1trial/(d1*d1) + n2trial/(d2*d2);
			double kyld = GetYield(mptr,np,delTime);
			double glam = 0.5*fnp12 - kyld*kyld/3.;
			double fnp1 = sqrt(fnp12);
			double slope = -(psKred*n1trial/(d1*d1*d1) + 2*Gred*n2trial/(d2*d2*d2)) - GetK2Prime(mptr,fnp1,delTime);
            
            // bisect if Newton out of range
            if( ((lambdak-xh)*slope-glam) * ((lambdak-xl)*slope-glam) >= 0. ||
               fabs(2.*glam) > fabs(dxold*slope) )
            {   dxold = dx;
                dx = 0.5*(xh-xl);
                lambdak = xl+dx;
                if(xl == lambdak) break;    // change in root is negligible
            }
            else
            {   dxold = dx;
                dx = glam/slope;
                double temp = lambdak;
                lambdak -= dx;
                if(temp == lambdak) break;  // change in root is negligible
            }
            
            // update and check convergence
            UpdateTrialAlpha(mptr,np,lambdak,(double)0.);
            if(LambdaConverged(step++,lambdak,dx)) break;
            
            // reset limits
            if(glam < 0.)
                xl = lambdak;
            else
                xh = lambdak;
        }
	}
	else
	{	while(true)
        {	// update iterative variables (lambda, alpha)
            double glam = strial - 2*Gred*lambdak - SQRT_TWOTHIRDS*GetYield(mptr,np,delTime);
            double slope = -2.*Gred - GetKPrime(mptr,np,delTime);
            
            // bisect if Newton out of range
            if( ((lambdak-xh)*slope-glam) * ((lambdak-xl)*slope-glam) >= 0. ||
                    fabs(2.*glam) > fabs(dxold*slope) )
            {   dxold = dx;
                dx = 0.5*(xh-xl);
                lambdak = xl+dx;
                if(xl == lambdak) break;    // change in root is negligible
            }
            else
            {   dxold = dx;
                dx = glam/slope;
                double temp = lambdak;
                lambdak -= dx;
                if(temp == lambdak) break;  // change in root is negligible
            }
            
            // update and check convergence
            UpdateTrialAlpha(mptr,np,lambdak,(double)0.);
            if(LambdaConverged(step++,lambdak,dx)) break;
            
            // reset limits
            if(glam < 0.)
                xl = lambdak;
            else
                xh = lambdak;
        }
	}
    
    // return final answer
    // cout << "   lambdak = " << (lambdak*SQRT_TWOTHIRDS/delTime) << endl;
	return lambdak;
}

/* Bracket the solution for lambda for safe Newton's method
    Subclass can override if have faster way to bracket
    ftrial is 3D or plane strain result for lamda=0 and it is positive
// Return lamNeg for f<0 (higher lambda) and lamPos where f>0 (lower lambda
*/
void IsoPlasticity::BracketSolution(MPMBase *mptr,int np,double strial,Tensor *stk,double delTime,
                                        double *lamNeg,double *lamPos)
{
    double epdot = 1.,gmax;
    int step=0;
    
    // take lambda = 0 as positive limit (to start)
    *lamPos = 0.;
    
    if(np==PLANE_STRESS_MPM)
	{	double n2trial = -stk->xx+stk->yy;
		n2trial *= 0.5*n2trial;
		n2trial += 2.*stk->xy*stk->xy;
		double n1trial = stk->xx+stk->yy;
		n1trial *= n1trial/6.;
        
        // find when plane stress term become negative
        while(step<20)
        {   // try above
            dalpha = epdot*delTime;
            alpint = mptr->GetHistoryDble() + dalpha;
            double lambdak = dalpha/SQRT_TWOTHIRDS;
			double d1 = (1 + psKred*lambdak);
			double d2 = (1.+2.*Gred*lambdak);
			double fnp12 = n1trial/(d1*d1) + n2trial/(d2*d2);
			double kyld = GetYield(mptr,np,delTime);
			gmax = 0.5*fnp12 - kyld*kyld/3.;
            if(gmax<0.) break;
            
            // update positive limit and go to next order of magnitude
            *lamPos = lambdak;
            epdot *= 10.;
            step++;
        }
    }
    else
    {   // find when strial 2 GRed sqrt(3/2) dalpha - sqrt(2/3)GetYield(alpha+dalpha,dalpha)
        // becomes negative
        while(step<20)
        {   // try above
            dalpha = epdot*delTime;
            alpint = mptr->GetHistoryDble() + dalpha;
            gmax = strial - 2*Gred*dalpha/SQRT_TWOTHIRDS - SQRT_TWOTHIRDS*GetYield(mptr,np,delTime) ;
            if(gmax<0.) break;
        
            // next block
            *lamPos = dalpha/SQRT_TWOTHIRDS;
            epdot *= 10.;
            step++;
        }
    }
    
    // exception if did not find answer in 20 orders of magnitude in strain rate
    if(step>=20)
        throw CommonException("Plasticity solution could not be bracketed","IsoPlasticity::BracketSolution");
    
    // set upper limits
    *lamNeg = dalpha/SQRT_TWOTHIRDS;
    
    //cout << "steps: " << step << ", epdot range: " << (*lamPos*SQRT_TWOTHIRDS/delTime) <<
    //        " to " << (*lamNeg*SQRT_TWOTHIRDS/delTime) << endl;
}

// decide if the numerical solution for lambda has converged
// subclass can override to change convergence rules
bool IsoPlasticity::LambdaConverged(int step,double lambda,double delLam)
{
	if(step>20 || fabs(delLam/lambda)<0.0001) return true;
	return false;
}

// material can override if history variable changes during elastic update (e.g., history dependent on plastic strain rate now zero)
void IsoPlasticity::ElasticUpdateFinished(MPMBase *mptr,int np,double delTime)
{
}

#pragma mark IsoPlaticity::Accessors

// archive material data for this material type when requested.
double IsoPlasticity::GetHistory(int num,char *historyPtr)
{
    double history=0.;
    if(num==1)
    {	double *cumStrain=(double *)historyPtr;
        history=*cumStrain;
    }
    return history;
}

// plastic strain needed to get deformation gradient for this material class
bool IsoPlasticity::HasPlasticStrainForGradient(void) { return TRUE; }

