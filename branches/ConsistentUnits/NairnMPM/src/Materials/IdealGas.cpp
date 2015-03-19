/********************************************************************************
    IdealGas.cpp
    nairn-mpm-fea

    Created by Edward Le and Peter Mackenzie on May 26, 2012.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/IdealGas.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "System/UnitsController.hpp"

#pragma mark IdealGas::Constructors and Destructors

// Constructors
IdealGas::IdealGas() {}

// Constructors with arguments 
IdealGas::IdealGas(char *matName) : HyperElastic(matName)
{
	P0   = -1.;			// required initial pressure
	rho  = -1.;			// required density (override default of 1)
	T0   = -1.;			// required initial temperature in Kelvin
}

#pragma mark IdealGas::Initialization

// Read material properties
char *IdealGas::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"P0")==0)
        return UnitsController::ScaledPtr((char *)&P0,gScaling,1.e6);
    
    else if(strcmp(xName,"T0")==0)
        return((char *)&T0);
	
    return(HyperElastic::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
const char *IdealGas::VerifyAndLoadProperties(int np)
{
	// make sure all were set
    if(P0 <= 0. || rho <= 0.0 || T0 <= 0.0 )
		return "Ideal gas material model needs positive parameters P0, rho, and T0";
	
	// Find ideal gas has Cv heat capacity in F-L/(mass-K)
	// For monotonic Ideal Gas, Cv = 1.5R for diatomic gas is 2.5R
	// If set to >1 is diatomic, otherwise monotonic (which is for not set too)
	if(heatCapacity>UnitsController::Scaling(1.e6))
	{	heatCapacity = 2.5*P0/(T0*rho);
		gammaAdiabatic = 7./5.;
	}
	else
	{	heatCapacity = 1.5*P0/(T0*rho);
		gammaAdiabatic = 5./3.;
	}
	CpMinusCv= P0/(T0*rho);
	
	// P0 in specific units (F/L^2 L^3/mass)
	P0sp=P0/rho;
	
    // call super class
    return HyperElastic::VerifyAndLoadProperties(np);
}

// if analysis not allowed, throw an exception
void IdealGas::ValidateForUse(int np) const
{	if(np==PLANE_STRESS_MPM)
	{	throw CommonException("IdealGas material cannot do 2D plane stress analysis",
							  "IdealGas::ValidateForUse");
	}
	
	if(thermal.reference<=0)
	{	throw CommonException("IdealGas material requires the simulation to set the stress free temperature in degrees K",
							  "IdealGas::ValidateForUse");
	}
	
	// call super class (why can't call super class?)
	HyperElastic::ValidateForUse(np);
}

// print mechanical properties output window
void IdealGas::PrintMechanicalProperties(void) const
{
	PrintProperty("P0",P0*UnitsController::Scaling(1.e-6),"");
	PrintProperty("T0",T0,"");
	cout << endl;
}

// If needed, a material can initialize particle state
// For example, ideal gas initializes to base line pressure
void IdealGas::SetInitialParticleState(MPMBase *mptr,int np) const
{
    // The initial state has particle mass mp = Vp * rho at T = thermal.reference
    // Imagine heating from T0 to T holding volume constant and find Kirchoff stress / rho0
    double Psp = -P0sp * (thermal.reference/T0);
    
    // set the particle pressure
	Tensor *sp=mptr->GetStressTensor();
	sp->xx = Psp;
	sp->yy = Psp;
	sp->zz = Psp;
    
    // Initial particle strains are zero (because J=1)
    
    // call super class for Cauchy Green strain
    HyperElastic::SetInitialParticleState(mptr,np);
}

#pragma mark IdealGas::Methods

// In unit nJ/(g-K)
double IdealGas::GetCpMinusCv(MPMBase *) const { return CpMinusCv; }

/* Take increments in strain and calculate new
    Particle: strains, rotation strain, stresses, strain energy, angle
    du are (gradient rates X time increment) to give deformation gradient change
*/
void IdealGas::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
    // Update strains and rotations and Left Cauchy strain
    // get determinent of incremental deformation gradient
    double detf = IncrementDeformation(mptr,du,NULL,np);
    
    // update stress
	Tensor *sp=mptr->GetStressTensor();
    double mPnsp = sp->xx;
    
	// compute specific pressure as p/rho = P0sp * (T/T0)
	// incrementally p(n+1)/rho = P0sp * (T(n+1)/T0)
	//							= P0sp * (Tn/T0) * (T(n+1)/Tn)
	//							= pn * (T(n+1)/Tn)
	// (note: pPreviousTemperature is T(n+1))
	double mPsp = (sp->xx/(1.-res->dT/mptr->pPreviousTemperature));

	// store in stress (which is minus the pressure)
	sp->xx = mPsp;
	sp->yy = mPsp;
	sp->zz = mPsp;
	
	// find the -P dV energy per unit mass dW/(rho0 V0) (nJ/g) as follows
    // dW/(rho0 V0) = - 0.5 * (pn+p(n+1))/rho0 * (V(n+1)-Vn)/V0, which simplifies to
    double dW = 0.5*(mPnsp*detf + mPsp)*(1.-1./detf);
    
    // this energy is tracked in work energy and no residual energy is tracked
    mptr->AddWorkEnergy(dW);
    
    // the same energy is tracked as heat (although it will be zero if adiabatic)
    // and is dissipated (which will cause heating if adiabatic
    // Update is Cv dT - dU
    IncrementHeatEnergy(mptr,res->dT,0.,dW);
}

#pragma mark IdealGas::Accessors

// return material type
const char *IdealGas::MaterialType(void) const { return "Ideal Gas (Hyperelastic)"; }

// Return the material tag
int IdealGas::MaterialTag(void) const { return IDEALGASMATERIAL; }

// Calculate current wave speed in L/sec.
double IdealGas::WaveSpeed(bool threeD,MPMBase *mptr) const
{	double Pspcurr;
	if(mptr!=NULL)
	{	Tensor *sp=mptr->GetStressTensor();
		Pspcurr = fmax(-sp->xx,0.);
	}
	else
		Pspcurr = P0sp;
	return sqrt(gammaAdiabatic*Pspcurr);
}



