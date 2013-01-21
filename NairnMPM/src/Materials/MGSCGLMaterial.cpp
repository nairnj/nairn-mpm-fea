/********************************************************************************
    MGSCGLMaterial.cpp
    NairnMPM
    
    Created by John Nairn, June 16, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
 
    This material is MGEOS only and can attach any hardening law. The name is
        code is form when it only allowed SCGLHardening law
********************************************************************************/

#include "Materials/MGSCGLMaterial.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Materials/HardeningLawBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"

#pragma mark MGSCGLMaterial::Constructors and Destructors

// Constructors
MGSCGLMaterial::MGSCGLMaterial() {}

// Constructors
MGSCGLMaterial::MGSCGLMaterial(char *matName) : IsoPlasticity(matName)
{
	gamma0=1.64;		// dimensionless
	C0=4004;			// m/sec
	S1=1.35;			// dimsionless
	S2=0.;				// dimsionless
	S3=0.;				// dimsionless
}

#pragma mark MGSCGLMaterial::Initialization

// Read material properties
char *MGSCGLMaterial::InputMat(char *xName,int &input)
{
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
	
    return(IsoPlasticity::InputMat(xName,input));
}

// verify settings and some initial calculations
const char *MGSCGLMaterial::VerifyProperties(int np)
{	
	// check properties (need here because IsotropicMat is skipped
	if(!read[G_PROP]) return "The shear modulus, G0, is missing";
	
	// needed because SetAnalysisProps never called
	CME3=betaI*concSaturation;

	// call plastic law, but skip IsotropicPlasticity and IsotropicMat
    const char *ptr = plasticLaw->VerifyProperties(np);
    if(ptr != NULL) return ptr;
	return MaterialBase::VerifyProperties(np);
}

// Constant properties used in constitutive law
void MGSCGLMaterial::InitialLoadMechProps(int makeSpecific,int np)
{	
    // hardening law properties
    plasticLaw->InitialLoadMechProps(makeSpecific,np);
	
    // Use in place of C0^2. Units are Pa cm^3/g such that get Pa when multiplied
    //      by a density in g/cm^3
    C0squared = 1000.*C0*C0;
	
    // Shear modulus with pressure dependence
	G0red = G*1.e6/rho;         // G0red = G/rho0
	Gred = G0red;               // Gred = G/rho = G rho0/(rho rho0) = J G0red
	Kred = C0squared;
    Keffred = C0squared;
	
	// thermal expansion is handled in EOS in GetPressureChange() so need to set
	// CTE3 used by Isoplasticity to 0;
	// Not sure if this material can handle solvent expansion?
	CTE3=0.;
	
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
    PrintProperty("G0",G,"");
	cout << endl;
    
	PrintProperty("S1",S1,"");
	PrintProperty("S2",S2,"");
	PrintProperty("S3",S3,"");
	cout << endl;
	
	// effective volumetric CTE (in ppm/K) alpha = rho0 gamma0 Cv / K
	double effAlpha = (1.e9*heatCapacityVol*gamma0)/C0squared;
	PrintProperty("a",effAlpha/3.,"");
	PrintProperty("T0",thermal.reference,"K");
	cout <<  endl;
    
    plasticLaw->PrintYieldProperties();

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
	{	throw CommonException("MGEOSMaterial material requires the simulation to set the stress free temperature in degrees K",
							  "MGSCGLMaterial::ValidateForUse");
	}
    
    if(np==PLANE_STRESS_MPM)
    {	throw CommonException("MGEOSMaterial material has not yet been updated to do plane stress calculations",
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
	
	// plane stress needs to adjust delV
	// Note: this calculation is based on initial state, because final state not known yet
	if(np==PLANE_STRESS_MPM)
    {   //Tensor *ep=mptr->GetStrainTensor();
        //double x0 = -(ep->xx+ep->yy+ep->zz);
		//Kred = k1 + x0*(2.*k2 + 3.*k3*x0);
		//double neta0=pow(1./(1.-x0),ONETHIRD);
		//Gred = G0red*(1.+GPpred*pressure0/neta0 + GTp*dTemp);
        Gred = G0red;
		delV *= 1./(Kred/(2.*Gred) + TWOTHIRDS);			// scale by (1-2nu)/(1-nu) for plane stress
	}
	
	// M-G EOS
    // Want specific pressure or pressure over current density (using J = rho0/rho)
    if(x>0.)
    {   // compression law
        // denominator = 1 - S1*x - S2*x^2 - S3*x^3
        double denom = 1./(1. - x*(S1 + x*(S2 + x*S3)));
    
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
    Gred = G0red * plasticLaw->GetShearRatio(mptr,pressure,J);
    
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

// This material is tracking specific Cauchy stress. On archiving need to know volume
// to get to actual Cauchy stress
double MGSCGLMaterial::GetCurrentRelativeVolume(MPMBase *mptr)
{   return mptr->GetRelativeVolume();
}

#pragma mark MGSCGLMaterial::Accessors

// Return the material tag
int MGSCGLMaterial::MaterialTag(void) { return MGEOSMATERIAL; }

// return unique, short name for this material
const char *MGSCGLMaterial::MaterialType(void) { return "MGEOS Material"; }

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


