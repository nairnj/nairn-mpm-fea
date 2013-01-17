/********************************************************************************
    MGSCGLMaterial.cpp
    NairnMPM
    
    Created by John Nairn, June 16, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/MGSCGLMaterial.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Global_Quantities/ThermalRamp.hpp"

#include "NairnMPM_Class/NairnMPM.hpp"

#pragma mark MGSCGLMaterial::Constructors and Destructors

// Constructors
MGSCGLMaterial::MGSCGLMaterial() {}

/* The default contructor should call a parent class constructor and
	then fill in any new initialization.
	*/
// Constructors
MGSCGLMaterial::MGSCGLMaterial(char *matName) : IsoPlasticity(matName)
{
	// defaults are some Tungsten properties
	yield=2200.;		// MPa
	yieldMax=0.;		// MPa
	beta=0.0;			// dimensionless
	nhard=1.0;			// dimensionless
	
	gamma0=1.64;		// dimensionless
	C0=4004;			// m/sec
	S1=1.35;			// dimsionless
	S2=0.;				// dimsionless
	S3=0.;				// dimsionless
	
	GPp=0.01e-3;		// MPa^-1
	GTp=-2.2e-4;		// K^-1
}

#pragma mark MGSCGLMaterial::Initialization

// Read material properties
char *MGSCGLMaterial::InputMat(char *xName,int &input)
{
	// unique properties here and
	//	yield found in Isoplasticity super class
	//	G found in IsotropicMat
	//	rho, Cv found in material base class
	
	// here are the rest
    if(strcmp(xName,"gamma0")==0)
    {	input=DOUBLE_NUM;
        return((char *)&gamma0);
    }
    
    else if(strcmp(xName,"C0")==0)
    {	input=DOUBLE_NUM;
        return((char *)&C0);
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
	
    else if(strcmp(xName,"GPpG0")==0)
    {	input=DOUBLE_NUM;
        return((char *)&GPp);
    }
	
    else if(strcmp(xName,"GTpG0")==0)
    {	input=DOUBLE_NUM;
        return((char *)&GTp);
    }
	
    else if(strcmp(xName,"betahard")==0)
    {	input=DOUBLE_NUM;
        return((char *)&beta);
    }
	
    else if(strcmp(xName,"nhard")==0)
    {	input=DOUBLE_NUM;
        return((char *)&nhard);
    }
	
    else if(strcmp(xName,"yieldMax")==0)
    {	input=DOUBLE_NUM;
        return((char *)&yieldMax);
    }
	
    return(IsoPlasticity::InputMat(xName,input));
}

// verify settings and some initial calculations
const char *MGSCGLMaterial::VerifyProperties(int np)
{	
	// check properties
	if(!readYield) return "The yield stress is missing";
	if(!read[G_PROP]) return "The shear modulus, G0, is missing";
	if(yieldMax<yield) return "The maximum yield stress is less than the initial yield stress";
	
	// needed because SetAnalysisProps never called
	CME3=betaI*concSaturation;

	// call super class
	return MaterialBase::VerifyProperties(np);
}

// Constant properties used in constitutive law
void MGSCGLMaterial::InitialLoadMechProps(int makeSpecific,int np)
{	
	yldred = yield*1.e6/rho;
	yldMaxred = yieldMax*1.e6/rho;
	
    // Use in place of C0^2. Units are Pa cm^3/g such that get Pa when multiplied
    //      by a density in g/cm^3
    C0squared = 1000.*C0*C0;
	
    // Shear modulus with pressure dependence
	G0red = G*1.e6/rho;         // G0red = G/rho0
	Gred = G0red;               // Gred = G/rho = G rho0/(rho rho0) = J G0red
	GPpred = GPp*rho*1.e-6;
	Kred = C0squared;
    Keffred = C0squared;
	
	// thermal expansion is handled in EOS in GetPressureChange() so need to set
	// CTE3 used by Isoplasticity to 0;
	// Not sure if this material can handle solvent expansion?
	CTE3=0.;
	
	// superclass sets yldred, but that replaced above, nothing else needed
	
	// plane stress correction done in GetPressureChange()
	psRed=1.;
}

// print mechanical properties to the results
void MGSCGLMaterial::PrintMechanicalProperties(void)
{
	// core properties
	PrintProperty("C0",C0,"m/s");
	PrintProperty("gam0",gamma0,"");
	PrintProperty("K",rho*Kred*1e-6,"");
	cout << endl;
    
	PrintProperty("S1",S1,"");
	PrintProperty("S2",S2,"");
	PrintProperty("S3",S3,"");
	cout << endl;
	
	// shear
	PrintProperty("G0",G,"");
	PrintProperty("Gp'/G0",GPp,"MPa^-1");
	PrintProperty("GT'/G0",GTp,"K^-1");
	cout << endl;
	
	// yield
	PrintYieldProperties();
	
	// effective CTE (in ppm/K) alpha = rho0 gamma0 Cv / K
	double effAlpha = (1.e9*heatCapacityVol*gamma0)/C0squared;
	PrintProperty("a",effAlpha/3.,"");
	PrintProperty("T0",thermal.reference,"K");
	cout <<  endl;
}

// print just yield properties to output window
void MGSCGLMaterial::PrintYieldProperties(void)
{
	PrintProperty("yld",yield,"");
	PrintProperty("beta",beta,"");
	PrintProperty("nhard",nhard,"");
	PrintProperty("yMax",yieldMax,"");
	cout << endl;
}

// Print transport properties
void MGSCGLMaterial::PrintTransportProperties(void)
{
	// Conductivity constants
	if(ConductionTask::active)
	{	MaterialBase::PrintTransportProperties();
		PrintProperty("Cv",heatCapacityVol,"J/(kg-K)");
	}
	else
	{	PrintProperty("Cp",heatCapacity,"J/(kg-K)");
		PrintProperty("Cv",heatCapacityVol,"J/(kg-K)");
	}
	cout << endl;
}

// if analysis not allowed, throw an exception
void MGSCGLMaterial::ValidateForUse(int np)
{	
	if(thermal.reference<=0)
	{	throw CommonException("MGSCGLMaterial material requires the simulation to set the stress free temperature in degrees K",
							  "MGSCGLMaterial::ValidateForUse");
	}
    
    if(np==PLANE_STRESS_MPM)
    {	throw CommonException("MGSCGLMaterial material has not yet been updated to do plane stress calculations",
                              "MGSCGLMaterial::ValidateForUse");
    }
	
	// call super class
	return IsoPlasticity::ValidateForUse(np);
}



#pragma mark MGSCGLMaterial::Custom Methods

/* To better handle large deformation in the M-G EOS, this method is overridden
    to calculate incremental deformation using large-deformation theory (i.e. from
    exponent of du), find delV, find Jnew, and then call back to Isoplasticity to
    finish shear parts using hypoelasticity
*/
void MGSCGLMaterial::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np)
{
    // Correct for swelling by finding total residual stretch (not incremental)
	double eres=0.,JresStretch=1.,deres=0.,dJresStretch=1.;
	if(DiffusionTask::active)
    {   eres += CME3*(mptr->pPreviousConcentration-DiffusionTask::reference);
        deres += CME3*DiffusionTask::dConcentration;
        JresStretch = (1.+eres)*(1.+eres)*(1.+eres);
        dJresStretch = (1.+deres)*(1.+deres)*(1.+deres);
    }
	
	// Get incremental deformation
	const Matrix3 dF = du.Exponential(incrementalDefGradTerms);

	// get new 2D deformation gradient to increment particle strain
    // plaine stress will need to add ep->zz when known
    Tensor *ep=mptr->GetStrainTensor();
    Tensor *eplast=mptr->GetPlasticStrainTensor();
    TensorAntisym *wrot = mptr->GetRotationStrainTensor();
    
    // get total strains
    double eptotxx = ep->xx + eplast->xx;
    double eptotxy = ep->xy + eplast->xy;
    double eptotyy = ep->yy + eplast->yy;
    double eptotzz = ep->zz + eplast->zz;
    
    if(np!=THREED_MPM)
    {   // previous particle deformation gradient
        const Matrix3 pF(1. + eptotxx,0.5*(eptotxy - wrot->xy),0.5*(eptotxy + wrot->xy),
                     1. + eptotyy,1. + eptotzz);
        
        // Trial deformation Gradient F = dF.pF
        const Matrix3 F = dF * pF;
        
        // redefine incremental deformation terms for large strain
        du = F - pF;
        
        // Trial update assuming elastic response
        double delV,Jnew;
        if(np==PLANE_STRESS_MPM)
        {   delV = psRed*(du(0,0)+du(1,1)-2.*eres)/3.;
            Jnew=1.;
        }
        else
        {   Jnew = F.determinant()/JresStretch;
            double dJ = dF.determinant()/dJresStretch;
            delV = (Jnew/3.)*(1.-1/dJ);
        }
        
        // Finish of shear parts and yield in the base IsoPlasticity class
        PlasticityConstLaw(mptr,du(0,0),du(1,1),du(0,1),du(1,0),du(2,2),delTime,np,delV,Jnew,eres);
    }
    else
    {   // two more strains
        double eptotxz = ep->xz + eplast->xz;
        double eptotyz = ep->yz + eplast->yz;
    
        // previous particle deformation gradient
        const Matrix3 pF(1. + eptotxx,             0.5*(eptotxy - wrot->xy), 0.5*(eptotxz - wrot->xz),
                     0.5*(eptotxy + wrot->xy), 1. + eptotyy,             0.5*(eptotyz - wrot->yz),
                     0.5*(eptotxz + wrot->xz), 0.5*(eptotyz + wrot->yz), 1. + eptotzz);
        
        // Trial deformation Gradient F = dF.pF
        const Matrix3 F = dF * pF;
        
        // redefine incremental deformation terms for large strain
        du = F - pF;
        
        // Trial update assuming elastic response
        double Jnew = F.determinant()/JresStretch;
        double dJ = dF.determinant()/dJresStretch;
        double delV = (Jnew/3.)*(1.-1/dJ);
        
        // Finish of shear parts and yield in the base IsoPlasticity class
        PlasticityConstLaw(mptr,du(0,0),du(1,1),du(2,2),du(0,1),du(1,0),du(0,2),du(2,0),
                           du(1,2),du(2,1),delTime,np,delV,Jnew,eres);
   }
}

// Find pressure using Mie-Gruneisen EOS
// Also get shear modulus here since it is not needed until after this is called
double MGSCGLMaterial::GetPressureChange(MPMBase *mptr,double &delV,double J,int np)
{
	// 3*delV is total incremental volumetric strain relative to free-swelling volume
    // J is total volume change - may need to reference to free-swelling volume if that works
	// Note that swelling looks like a problem because the sums of strains needs to be adjusted
	//		to stress-free state
	
	// compression
    double x = 1.-J;        // current compression J = 1-x
	Tensor *sptr=mptr->GetStressTensor();
	double pressure0 = -(sptr->xx+sptr->yy+sptr->zz)/3.;
	double dTemp=mptr->pPreviousTemperature-thermal.reference;
	
	// plane stress needs to adjust delV
	// Note: this calculation is based on initial state, because final state not known yet
	if(np==PLANE_STRESS_MPM)
    {   Tensor *ep=mptr->GetStrainTensor();
        double x0 = -(ep->xx+ep->yy+ep->zz);
		//Kred = k1 + x0*(2.*k2 + 3.*k3*x0);
		double neta0=pow(1./(1.-x0),ONETHIRD);
		Gred = G0red*(1.+GPpred*pressure0/neta0 + GTp*dTemp);
		delV *= 1./(Kred/(2.*Gred) + TWOTHIRDS);			// scale by (1-2nu)/(1-nu) for plane stress
	}
	
	// M-G EOS
    // Want specific pressure or pressure over current density (using J = rho0/rho)
    if(x>0.)
    {   // compression law
        // denominator = 1 - S1*x - S2*x^2 - S2^x^3
        double denom = 1./(1. - x*(S1 +x*(S2 + x*S3)));
    
        // current effective and reduced (by rho0) bulk modulus
        Keffred = C0squared*(1.-0.5*gamma0*x)*denom*denom;
    }
    else
    {   // In tension use low-strain bulk modulus
        Keffred = C0squared;
    }
    
    // Pressure from bulk modulus and an energy term
    double e = mptr->GetStrainPlusPlastEnergy()
                + 1000.*heatCapacityVol*(mptr->pPreviousTemperature - thermal.reference*exp(gamma0*x));
	double pressure = J*(Keffred*x + gamma0*e);
    
	// SCGL shear modulus and save Gratio = J G/G0 for later calculations
    // Note: J in Gred and Gratio is so that where they are used, they give
    //          specific Cauchy stress
	double neta=pow(1./J,ONETHIRD);
    Gratio = J*(1.+GPpred*pressure/neta + GTp*dTemp);
	Gred = G0red*Gratio;
    
	// plane stress terms
	if(np==PLANE_STRESS_MPM)
	{	//Kred = k1 + x*(2.*k2 + 3.*k3*x);
		double psRed0 = 1./(Kred/(2.*Gred) + TWOTHIRDS);	// (1-2nu)/(1-nu) for plane stress
		psLr2G = (Kred/(2.*Gred) - ONETHIRD)*psRed0;		// nu/(1-nu) to find ezz
		psKred = Kred*psRed0;                               // E/(3(1-v)) to find lambda
	}
	
	// return change in pressure
	return pressure-pressure0;
}

// Return yield stress for current conditions and it is specific Cauchy stress
//   (alpint for cum. plastic strain and dalpha/delTime for plastic strain rate)
// yield = yldred*(1 + beta ep)^n * Gred/G0red, where ep=alpint
// but leading term is limited to yldMaxred
double MGSCGLMaterial::GetYield(MPMBase *mptr,int np,double delTime)
{	return fmin(yldred*pow(1.+beta*alpint,nhard),yldMaxred)*Gratio;
}

// Get derivative of sqrt(2./3.)*yield with respect to lambda for plane strain and 3D
// ... and using dep/dlambda = sqrt(2./3.)
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)lamda or depdot/dlambda = sqrt(2./3.)/delTime
// ... and as specfic Cauchy stress
double MGSCGLMaterial::GetKPrime(MPMBase *mptr,int np,double delTime)
{	
	// slope zero if in constant max yield condition
	if(yldred*pow(1.+beta*alpint,nhard)>=yldMaxred) return 0.;
	
	// return slope
	double factor=yldred*Gratio;
	double bfactor = DbleEqual(nhard,1.) ? beta :
				beta*nhard*pow(1.+beta*alpint,nhard-1.) ;
	return TWOTHIRDS*factor*bfactor;
}

// this material does not support plane stress calculations
double MGSCGLMaterial::GetK2Prime(MPMBase *mptr,double fnp1,double delTime)
{
	// slope zero if in constant max yield condition
	if(yldred*pow(1.+beta*alpint,nhard)>=yldMaxred) return 0.;

	double factor=yldred*Gratio;
	return SQRT_EIGHT27THS*factor*factor*beta*nhard*pow(1.+beta*alpint,2.*nhard-1)*fnp1;
}

// This material is tracking specific Cauchy stress. On archiving need to know volume
// to get to actual Cauchy stress
double MGSCGLMaterial::GetCurrentRelativeVolume(MPMBase *mptr)
{   return mptr->GetRelativeVolume();
}

#pragma mark MGSCGLMaterial::Accessors

// Return the material tag
int MGSCGLMaterial::MaterialTag(void) { return MGSCGLMATERIAL; }

// return unique, short name for this material
const char *MGSCGLMaterial::MaterialType(void) { return "M-G, SGCL Material"; }

/*	calculate wave speed in mm/sec. Uses initial sqrt((K+4G/3)/rho) which is dilational wave speed
    K in Pa is 1000*rho*C0^2, G is in MPa, rho is in g/cm^3
*/
double MGSCGLMaterial::WaveSpeed(bool threeD,MPMBase *mptr) { return 1000.*sqrt(C0*C0+4000.*G/(3.*rho)); }

/*	calculate current wave speed in mm/sec. Uses sqrt((K+4G/3)/rho) which is dilational wave speed
    but K and G are current values of rho0*Keffred and rho0*Gred/J (in Pa) so were want
    1000 sqrt(rho0(Keffred+4Gred/(3J))/(1000 rho0/J))
 */
double MGSCGLMaterial::CurrentWaveSpeed(bool threeD,MPMBase *mptr)
{
    double J = mptr->GetRelativeVolume();
    return 1000.*sqrt((J*Keffred + 4.*Gred/3.)/1000.);
}


