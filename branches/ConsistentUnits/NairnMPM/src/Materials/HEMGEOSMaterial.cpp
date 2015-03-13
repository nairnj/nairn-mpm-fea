/********************************************************************************
	MGSCGLMaterial.cpp
	nairn-mpm-fea

	Created by John Nairn, Feb 18, 2013.
	Copyright (c) 2013 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/HEMGEOSMaterial.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Materials/HardeningLawBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"

#pragma mark HEMGEOSMaterial::Constructors and Destructors

// Constructors
HEMGEOSMaterial::HEMGEOSMaterial() {}

// Constructors
HEMGEOSMaterial::HEMGEOSMaterial(char *matName) : HEIsotropic(matName)
{
	gamma0=1.64;		// dimensionless
	C0=4004;			// m/sec
	S1=1.35;			// dimsionless
	S2=0.;				// dimsionless
	S3=0.;				// dimsionless
	Kbulk = 1.;			// not used
}

#pragma mark HEMGEOSMaterial::Initialization

// Read material properties
char *HEMGEOSMaterial::InputMaterialProperty(char *xName,int &input,double &gScaling)
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
	
   return(HEIsotropic::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
const char *HEMGEOSMaterial::VerifyAndLoadProperties(int np)
{
	// call plastic law first
    const char *ptr = plasticLaw->VerifyAndLoadProperties(np);
    if(ptr != NULL) return ptr;
    
	// check properties (need here because IsotropicMat is skipped
	if(G1<0) return "The shear modulus, G1, is missing";
	
	// MU in specific units using initial rho
	// for MPM (units N/m^2 cm^3/g)
	G1sp = G1*1.0e+06/rho;
	
    // Use in place of C0^2. Units are Pa cm^3/g such that get Pa when multiplied
    //      by a density in g/cm^3
	// Equal to reduced bulk modulus
    Ksp = C0squared = 1000.*C0*C0;
	
    // Shear modulus with pressure dependence
	Kbulk = 1e-6*rho*C0squared;			// initial bulk modulus in MPa to print
	
	// expansion coefficients - affect on pressure is handled by eos, but
	// needed for shear parts (which in large deformation have Jres component)
	double effAlpha = (1.e9*heatCapacity*gamma0)/C0squared;
	CTE1 = 1.e-6*effAlpha/3.;
	
	// this material not coupled to moisture expansion
	betaI = 0.;
	CME1 = 0.;;
    
    // for Cp-Cv
    Ka2sp = 0.001*Ksp*CTE1*CTE1;
	
	// skip Hyperelstic methods
	return MaterialBase::VerifyAndLoadProperties(np);
}

// print mechanical properties to the results
void HEMGEOSMaterial::PrintMechanicalProperties(void) const
{
	// core properties
	PrintProperty("C0",C0,"m/s");
	PrintProperty("gam0",gamma0,"");
	PrintProperty("K",Kbulk,"");
    PrintProperty("G1",G1,"");
	cout << endl;
    
	PrintProperty("S1",S1,"");
	PrintProperty("S2",S2,"");
	PrintProperty("S3",S3,"");
	cout << endl;
	
	// effective volumetric CTE (in ppm/K) alpha = rho0 gamma0 Cv / K
	double effAlpha = (1.e9*heatCapacity*gamma0)/C0squared;
	PrintProperty("a",effAlpha/3.,"");
	PrintProperty("T0",thermal.reference,"K");
	cout <<  endl;
    
    plasticLaw->PrintYieldProperties();
	
	// skip super class, but call it's superclass
	HyperElastic::PrintMechanicalProperties();
}

// Print transport properties
void HEMGEOSMaterial::PrintTransportProperties(void) const
{
	// Conductivity constants
	if(ConductionTask::active)
	{	MaterialBase::PrintTransportProperties();
	}
	else if(!ConductionTask::adiabatic)
	{	PrintProperty("Cv",heatCapacity,"J/(kg-K)");
		cout << endl;
	}
}

// if analysis not allowed, throw an exception
void HEMGEOSMaterial::ValidateForUse(int np) const
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
	HEIsotropic::ValidateForUse(np);
}


#pragma mark MGSCGLMaterial::Custom Methods

// Get plastic properties or NULL on memory error
void *HEMGEOSMaterial::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer) const
{
	HEPlasticProperties *p = (HEPlasticProperties *)matBuffer;
 	p->hardProps = plasticLaw->GetCopyOfHardeningProps(mptr,np,altBuffer);
	// Gred and Kred found in UpdatePressure() - do not use before that
	return p;
}

// This method handles the pressure equation of state. Its tasks are
// 1. Calculate the new pressure
// 2. Update particle pressure
// 3. Increment the particle energy due to dilation
// 4. Call plasticLaw to see if it wants to change the shear modulus
// J = V(T,c)/V0(T,c), Jeff = V(T,c)/V0(Tref,cref), Jres = V0(T,c)/V0(Tref,cref) = J/Jeff
// Jtot = V(T,c)/V0(Trec,cref), Jres = V0(T,c)/V0(Tref,cref) for free expansion, J = V(T,c)/V0(T,c)
// Jn+1 = (detdF/detdFres) Jn, Jresn+1 = detdFres Jresn, Jtot = detdF Jtotn
// detdFres = (1+dres)^3 (approximately)
// Here Tref and cref are starting conditions and T and c are current temperature and moisture
void HEMGEOSMaterial::UpdatePressure(MPMBase *mptr,double J,double detdF,int np,double Jeff,
									 double delTime,HEPlasticProperties *p,ResidualStrains *res,double detdFres) const
{
    // J is total volume change - may need to reference to free-swelling volume if that works
	// Note that swelling looks like a problem because the sums of strains needs to be adjusted
	//		to stress-free state
	
	// previous pressure
	double P,P0 = mptr->GetPressure();
	double delV = 1. - 1./detdF;
	double dTq0;
    
    // M-G EOS
    // Want specific pressure or pressure over current density (using J = rho0/rho)
    if(J>1.)
	{	// new compression J(k+1) = 1-x(k+1)
		double x = 1.-J;
		
		// compression law
        // denominator = 1 - S1*x - S2*x^2 - S3*x^3
        double denom = 1./(1. - x*(S1 + x*(S2 + x*S3)));
		
        // law not valid if denominator passes zero
        if(denom<0)
        {   cout << "# Excessive x = " << x << endl;
            mptr->Describe();
        }
        
        // current effective and reduced (by rho0) bulk modulus
        p->Kred = C0squared*(1.-0.5*gamma0*x)*denom*denom;
		
		// Pressure from bulk modulus and an energy term
		double e = mptr->GetInternalEnergy();
		P = J*(p->Kred*x + gamma0*e);
		
		// particle isentropic temperature increment
		dTq0 = -J*gamma0*mptr->pPreviousTemperature*delV;
    }
    else
    {   // In tension hyperelastic law P = - K0(J-1)
        p->Kred = C0squared*Jeff;
		P = -J*C0squared*(Jeff-1.);
		//P = P0 - J*p->Kred*delV;
		
		// particle isentropic temperature increment
		double Kratio = Jeff;
		dTq0 = -J*Kratio*gamma0*mptr->pPreviousTemperature*delV;
    }
	
    // artifical viscosity
	// delV is total incremental volumetric strain = total Delta(V)/V
    double QAVred = 0.,AVEnergy=0.;
    if(delV<0. && artificialViscosity)
    {   double c = sqrt(p->Kred*J/1000.);        // m/sec
        QAVred = GetArtificalViscosity(delV/delTime,c);
        if(ConductionTask::AVHeating) AVEnergy = fabs(QAVred*delV);
    }
    
    // set final pressure
    mptr->SetPressure(P+QAVred);
    
    // work energy is dU = -P dV + s.de(total)
	// Here do hydrostatic terms, deviatoric later
    double avgP = 0.5*(P0+P);
	double delVres = 1. - 1./detdFres;
    mptr->AddWorkEnergyAndResidualEnergy(-avgP*delV,-avgP*delVres);
    
    // heat energy is Cv (dT - dTq0) - dPhi - |QAVred*delV|
	// Here do Cv (dT - dTq0) - |QAVred*delV| term and dPhi is done later
    IncrementHeatEnergy(mptr,res->dT,dTq0,AVEnergy);
	
	// SCGL and SL shear modulus and save Gratio = Jeff G/G0 for later calculations
    // Note: Jeff in Gred and Gratio is so that where they are used, they give
    //          specific Cauchy stress
    p->Gred = G1sp * plasticLaw->GetShearRatio(mptr,avgP,Jeff,p->hardProps);
}

#pragma mark MGSCGLMaterial::Accessors

// Return the material tag
int HEMGEOSMaterial::MaterialTag(void) const { return HEMGEOSMATERIAL; }

// return unique, short name for this material
const char *HEMGEOSMaterial::MaterialType(void) const { return "Hyperelastic MGEOS Material"; }

/*	calculate current wave speed in mm/sec. Uses sqrt((K+4G/3)/rho) which is dilational wave speed
 but K and G are current values of rho0*Keffred and rho0*Gred/J (in Pa) so were want
 1000 sqrt(rho0(Keffred+4Gred/(3J))/(1000 rho0/J))
 */
double HEMGEOSMaterial::CurrentWaveSpeed(bool threeD,MPMBase *mptr) const
{
    // compressive volumetric strain x = 1-J
    double J = mptr->GetRelativeVolume();
    double x = 1. - J;
    
    // get K/rho0, but this ignores slope of energy term
    double KcurrRed;
    if(x>0.)
    {   // compression law
        // denominator = 1 - S1*x - S2*x^2 - S3*x^3
        double denom = 1./(1. - x*(S1 + x*(S2 + x*S3)));
        
        // current effective and reduced (by rho0) bulk modulus
        KcurrRed = C0squared*(1.-0.5*gamma0*x)*denom*denom;
    }
    else
    {   // In tension use low-strain bulk modulus
        KcurrRed = C0squared;
    }
    KcurrRed *= J;          // convert to K/rho
    
    // get G/rho at current pressure
    double e = mptr->GetInternalEnergy();
    double pressure = J*(KcurrRed*x + gamma0*e);
	// MUST FIX
    double GcurrRed = G1sp * plasticLaw->GetShearRatio(mptr,pressure,J,NULL);
    
    // return current save speed
    return 1000.*sqrt((KcurrRed + 4.*GcurrRed/3.)/1000.);
}


