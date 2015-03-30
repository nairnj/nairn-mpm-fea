/********************************************************************************
    MGSCGLMaterial.cpp
    nairn-mpm-fea
    
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
#include "System/UnitsController.hpp"

#pragma mark MGSCGLMaterial::Constructors and Destructors

// Constructors
MGSCGLMaterial::MGSCGLMaterial() {}

// Constructors
MGSCGLMaterial::MGSCGLMaterial(char *matName) : IsoPlasticity(matName)
{
	gamma0=1.64;		// dimensionless
	C0=4004000.;		// mm/sec
	S1=1.35;			// dimsionless
	S2=0.;				// dimsionless
	S3=0.;				// dimsionless
}

#pragma mark MGSCGLMaterial::Initialization

// Read material properties
char *MGSCGLMaterial::InputMaterialProperty(char *xName,int &input,double &gScaling)
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
	
    return(IsoPlasticity::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
const char *MGSCGLMaterial::VerifyAndLoadProperties(int np)
{	
	// check properties (need here because IsotropicMat is skipped
	if(!read[G_PROP]) return "The shear modulus, G0, is missing";
	
	// needed because SetAnalysisProps never called
	CME3=betaI*concSaturation;

	// call plastic law, but skip IsotropicPlasticity and IsotropicMat
    const char *ptr = plasticLaw->VerifyAndLoadProperties(np);
    if(ptr != NULL) return ptr;
	
    // Use in place of C0^2. Units are L^2/sec^2 = F/L^2 L^3/mass
	// Equal to reduced bulk modulus
    C0squared = C0*C0;
 	
    // Shear modulus with pressure dependence
	G0red = G/rho;					// G0red = G/rho0
	pr.Gred = G0red;				// Gred = G/rho = G rho0/(rho rho0) = J G0red
	pr.Kred = C0squared;
	
	// thermal expansion is handled in EOS in UpdatePressure() so need to set
	// CTE3 used by Isoplasticity to 0;
	// Not sure if this material can handle solvent expansion?
	CTE3=0.;
	
	// skip to base class
	return MaterialBase::VerifyAndLoadProperties(np);
}

// print mechanical properties to the results
void MGSCGLMaterial::PrintMechanicalProperties(void) const
{
	// core properties
	PrintProperty("C0",C0*UnitsController::Scaling(1.e-3),UnitsController::Label(ALTVELOCITY_UNITS));
	PrintProperty("gam0",gamma0,"");
	PrintProperty("K",rho*pr.Kred*UnitsController::Scaling(1.e-6),"");
    PrintProperty("G0",G*UnitsController::Scaling(1.e-6),"");
	cout << endl;
    
	PrintProperty("S1",S1,"");
	PrintProperty("S2",S2,"");
	PrintProperty("S3",S3,"");
	cout << endl;
	
	// effective volumetric CTE (in ppm/K) alpha = rho0 gamma0 Cv / K
	double effAlpha = 1.e6*(heatCapacity*gamma0/C0squared);
	PrintProperty("a",effAlpha/3.,"");
	PrintProperty("T0",thermal.reference,"K");
	cout <<  endl;
    
    plasticLaw->PrintYieldProperties();
}

// Print transport properties
void MGSCGLMaterial::PrintTransportProperties(void) const
{
	// Conductivity constants
	if(ConductionTask::active)
	{	MaterialBase::PrintTransportProperties();
	}
	else if(!ConductionTask::adiabatic)
	{	PrintProperty("C",heatCapacity*UnitsController::Scaling(1.e-6),UnitsController::Label(HEATCAPACITY_UNITS));
		cout << endl;
	}
}

// if analysis not allowed, throw an exception
void MGSCGLMaterial::ValidateForUse(int np) const
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
	IsoPlasticity::ValidateForUse(np);
}


#pragma mark MGSCGLMaterial::Custom Methods

// Isotropic material can use read-only initial properties
void *MGSCGLMaterial::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer) const
{
	PlasticProperties *p = (PlasticProperties *)matBuffer;
	*p = pr;
 	p->hardProps = plasticLaw->GetCopyOfHardeningProps(mptr,np,altBuffer);
	
	// Gratio and other properties loaded later
	return p;
}

/* To better handle large deformation in the M-G EOS, this method is overridden
    to calculate incremental deformation using large-deformation theory (i.e. from
    exponent of du), find delV, find Jnew, and then call back to Isoplasticity to
    finish shear parts using hypoelasticity
*/
void MGSCGLMaterial::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
    // Correct for swelling by finding total residual stretch (not incremental)
	double eresTot=0.,JresStretch=1.,eres=0.,detdFres=1.;
	if(DiffusionTask::active)
    {   eresTot += CME3*(mptr->pPreviousConcentration-DiffusionTask::reference);
        eres += CME3*res->dC;
        JresStretch = (1.+eresTot)*(1.+eresTot)*(1.+eresTot);
		detdFres = (1.+eres)*(1.+eres)*(1.+eres);
    }
	
	// Get incremental deformation
	const Matrix3 dF = du.Exponential(incrementalDefGradTerms);
    
    // get previous deformation
    const Matrix3 pF = mptr->GetDeformationGradientMatrix();

    // Trial deformation Gradient F = dF.pF
    const Matrix3 F = dF * pF;
    
    // reassign incremental strains
    du = F - pF;
    
    double detdF = dF.determinant();                    // = V(k+1)/V(k)
    double delVtot = 1. - 1./detdF;                     // = (V(k+1)-V(k))/V(k+1)
    double delV = 1. - detdFres/detdF;                  // = (V(k+1)/Vsf(k+1) - V(k)/Vsf(k))/(V(k+1)/Vsf(k+1))
    double Jnew = F.determinant()/JresStretch;          // = V(k+1)/Vsf(k+1)
	PlasticProperties *p = (PlasticProperties *)properties;
    p->delVLowStrain = du.trace() -  3.*eres;
	
	// SCGL and SL chnage shear modulus, here Gratio =  J G/G0
    // Note: J in Gred and Gratio is so that where they are used, they give
    //          specific Cauchy stress
	double Gratio = plasticLaw->GetShearRatio(mptr,mptr->GetPressure(),Jnew,p->hardProps);
	p->Gred = G0red*Gratio;
    
    // artifical viscosity
    p->QAVred = 0.;
    if(delVtot<0. && artificialViscosity)
	{	p->QAVred = GetArtificalViscosity(delVtot/delTime,sqrt(p->Kred*Jnew*JresStretch));
    }
    
    if(np!=THREED_MPM)
    {   // Finish of shear parts and yield in the base IsoPlasticity class
        PlasticityConstLaw(mptr,du(0,0),du(1,1),du(0,1),du(1,0),du(2,2),delTime,np,delV,Jnew,eres,p,res);
    }
    else
    {   // Finish of shear parts and yield in the base IsoPlasticity class
        PlasticityConstLaw(mptr,du(0,0),du(1,1),du(2,2),du(0,1),du(1,0),du(0,2),du(2,0),
                           du(1,2),du(2,1),delTime,np,delV,Jnew,eres,p,res);
    }
}

// This method handles the pressure equation of state. Its tasks are
// 1. Calculate the new pressure
// 2. Update particle pressure
// 3. Increment the particle energy due to dilation
// 4. Call plasticLaw to see if it wants to change the shear modulus
// 5. Change delV to low strain result for subsequent plasticity calcs
// Notes:
//  delV is relative incremental volume change on this step = (V(k+1)-V(k))/V(k+1)
//  J is total volume change at end of step
void MGSCGLMaterial::UpdatePressure(MPMBase *mptr,double &delV,double J,int np,PlasticProperties *p,ResidualStrains *res,double eres) const
{
	// delV is total incremental volumetric strain relative to free-swelling volume
    // J is total volume change - may need to reference to free-swelling volume if that works
	// Note that swelling looks like a problem because the sums of strains needs to be adjusted
	//		to stress-free state
	
	// compression
    double x = 1.-J;        // new compression J(k+1) = 1-x(k+1)
    
    // M-G EOS
    // Want specific pressure or pressure over current density (using J = rho0/rho)
    if(x>0.)
    {   // compression law
        // denominator = 1 - S1*x - S2*x^2 - S3*x^3
        double denom = 1./(1. - x*(S1 + x*(S2 + x*S3)));

        // law not valid if denominator passes zero
        if(denom<0)
        {   cout << "#Excessive x = " << x << endl;
            mptr->Describe();
        }
        
        // current effective and reduced (by rho0) bulk modulus
        p->Kred = C0squared*(1.-0.5*gamma0*x)*denom*denom;
    }
    else
    {   // In tension use P = -K(J-1)
        p->Kred = C0squared;
    }
	
	// update plane stress terms now
	if(np==PLANE_STRESS_MPM)
	{	// these are terms for plane stress calculations only
		p->psRed = 1./(p->Kred/(2.*p->Gred) + 2./3.);					// (1-2nu)/(1-nu) for plane stress
		p->psLr2G = (p->Kred/(2.*p->Gred) - 1./3.)*p->psRed;			// nu/(1-nu) to find ezz
		p->psKred = p->Kred*p->psRed;									// E/(3(1-v)) to find lambda
	}

    // Pressure from bulk modulus and an energy term
    double e = mptr->GetInternalEnergy();
	double P0 = mptr->GetPressure();
    double P = J*(p->Kred*x + gamma0*e);
    
    // set final pressure
    mptr->SetPressure(P+p->QAVred);
    
    // particle isentropic temperature increment
    delV += 3.*eres;        // need total volume change now (better to use delVtot, but not available here without more calcs)
    double dTq0 = -gamma0*mptr->pPreviousTemperature*J*delV;
    
    // work energy is dU = -P dV + s.de(total)
	// Here do hydrostatic terms, deviatoric later
    double avgP = 0.5*(P0+P);
    mptr->AddWorkEnergyAndResidualEnergy(-avgP*delV,-3.*avgP*eres);
    
    // heat energy is Cv (dT - dTq0) - dPhi - QAVred*delV
	// Here do Cv (dT - dTq0) - QAVred*delV term and dPhi is done later
    double AVEnergy = ConductionTask::AVHeating ? fabs(p->QAVred*delV) : 0. ;
    IncrementHeatEnergy(mptr,res->dT,dTq0,AVEnergy);
    
    // reset for small-strain plasticity
    delV = p->delVLowStrain;
    
}

// This material is tracking specific Cauchy stress. On archiving need to know volume
// to get to actual Cauchy stress
double MGSCGLMaterial::GetCurrentRelativeVolume(MPMBase *mptr) const
{  
    return mptr->GetRelativeVolume();
}

#pragma mark MGSCGLMaterial::Accessors

// convert J to K using isotropic method
Vector MGSCGLMaterial::ConvertJToK(Vector d,Vector C,Vector J0,int np)
{	double KLS = rho*pr.Kred;
	double nuLS = (3.*KLS-2.*G)/(6.*KLS+2.*G);
	return IsotropicJToK(d,C,J0,np,nuLS,G);
}

// Return the material tag
int MGSCGLMaterial::MaterialTag(void) const { return MGEOSMATERIAL; }

// return unique, short name for this material
const char *MGSCGLMaterial::MaterialType(void) const { return "MGEOS Material"; }

// calculate wave speed in mm/sec. Uses initial sqrt((K+4G/3)/rho) which is dilational wave speed
// K/rho is C0squared in mm^2/sec^2, G is in Pa, rho is in g/mm^3
double MGSCGLMaterial::WaveSpeed(bool threeD,MPMBase *mptr) const { return sqrt(C0squared+4.*G/(3.*rho)); }

// Calculate current wave speed in mm/sec. Uses sqrt((K+4G/3)/rho) which is dilational wave speed
// but K/rho = Kred*J and G/rho = Gred*J (in mm^2/sec^2)
double MGSCGLMaterial::CurrentWaveSpeed(bool threeD,MPMBase *mptr) const
{
    // compressive volumetric strain x = 1-J
    double J = mptr->GetRelativeVolume();
    double x = 1. - J;
    
    // get K/rho0, but this ignores slope of energy term
    double KcurrRed,pressure;
    if(x>0.)
    {   // compression law
        // denominator = 1 - S1*x - S2*x^2 - S3*x^3
        double denom = 1./(1. - x*(S1 + x*(S2 + x*S3)));
        
        // current effective and reduced (by rho0) bulk modulus
        KcurrRed = C0squared*(1.-0.5*gamma0*x)*denom*denom;
		double e = mptr->GetInternalEnergy();
		pressure = J*(KcurrRed*x + gamma0*e);
    }
    else
    {   // In tension use low-strain bulk modulus
        KcurrRed = C0squared;
		pressure = -J*KcurrRed*(J-1.);
    }
    KcurrRed *= J;          // convert to K/rho
    
    // get G/rho at current pressure
    double GcurrRed = J*(G0red * plasticLaw->GetShearRatio(mptr,pressure,J,NULL));
    
    // return current save speed
    return sqrt((KcurrRed + 4.*GcurrRed/3.));
}

// if a subclass material supports artificial viscosity, override this and return TRUE
bool MGSCGLMaterial::SupportsArtificialViscosity(void) const { return true; }

