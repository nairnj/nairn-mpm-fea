/********************************************************************************
 IdealGas.cpp
 NairnMPM
 
 Created by Edward Le and Peter Mackenzie on May 26, 2012.
 Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/IdealGas.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Exceptions/CommonException.hpp"
 
#pragma mark IdealGas::Constructors and Destructors

// Constructors
IdealGas::IdealGas() {}

// Constructors with arguments 
IdealGas::IdealGas(char *matName) : HyperElastic(matName)
{
	P0   = -1.;			// required initial pressure in MPa
	rho  = -1.;			// required density (override default of 1) in g/cm^3
	T0   = -1.;			// required initial temperature in Kelvin
}

#pragma mark IdealGas::Initialization

// print mechanical properties output window
void IdealGas::PrintMechanicalProperties(void)
{
	PrintProperty("P0", P0, "");
	PrintProperty("T0", T0, "");
	cout << endl;
}
	
// Read material properties
char *IdealGas::InputMat(char *xName,int &input)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"P0")==0)
        return((char *)&P0);
    
    else if(strcmp(xName,"T0")==0)
        return((char *)&T0);
	
    return(HyperElastic::InputMat(xName,input));
}

// verify settings and some initial calculations
const char *IdealGas::VerifyProperties(int np)
{
	// make sure all were set
    if(P0 <= 0. || rho <= 0.0 || T0 <= 0.0 )
		return "Ideal gas material model needs positive parameters P0, rho, and T0";
	
	// Find ideal gas has heat capacity in J/(kg-K) = mJ/(g-K)
	// Must set C to Cv to work well with theory, even though output and equations
	//    for conductivity claim Cp is needed. For monotonic Ideal Gas, Cv = 1.5R
	//	  for diatomic is is 2.5R
	// If set to >1 is diatomic, otherwise monotonic (which is for not set too)
	if(heatCapacity>1.)
		heatCapacity = 2500.*P0/(T0*rho);
	else
		heatCapacity = 1500.*P0/(T0*rho);
	
    // call super class
    return HyperElastic::VerifyProperties(np);
}

// if analysis not allowed, throw an exception
void IdealGas::ValidateForUse(int np)
{	if(np==PLANE_STRESS_MPM)
	{	throw CommonException("IdealGas material cannot do 2D plane stress analysis",
							  "IdealGas::ValidateForUse");
	}
	
	if(thermal.reference<=0)
	{	throw CommonException("IdealGas material requires the simulation to set the stress free temperature in degrees K",
							  "IdealGas::ValidateForUse");
	}
	
	// call super class (why can't call super class?)
	return HyperElastic::ValidateForUse(np);
}

// Private properties used in constitutive law
void IdealGas::InitialLoadMechProps(int makeSpecific,int np)
{
	// P0 in specific units for MPM of N/m^2 cm^3/g
	P0sp=P0*1.0e+06/rho;
}

// If needed, a material can initialize particle state
// For example, ideal gas initializes to base line pressure
void IdealGas::SetInitialParticleState(MPMBase *mptr,int np)
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
    return HyperElastic::SetInitialParticleState(mptr,np);
}


#pragma mark IdealGas::Methods

/* Take increments in strain and calculate new
    Particle: strains, rotation strain, stresses, strain energy, angle
    du are (gradient rates X time increment) to give deformation gradient change
*/
void IdealGas::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np)
{
    // Update strains and rotations and Left Cauchy strain
    // get determinent of incremental deformation gradient
    double detf = IncrementDeformation(mptr,du,NULL,np);;
    
    // update stress
	Tensor *sp=mptr->GetStressTensor();
    double mPnsp = sp->xx;
    
	// compute specific pressure as p/rho = P0sp * (T/T0)
	// incrementally p(n+1)/rho = P0sp * (T(n+1)/T0)
	//							= P0sp * (Tn/T0) * (T(n+1)/Tn)
	//							= pn * (T(n+1)/Tn)
	// (note: pPreviousTemperature is T(n+1))
	double mPsp = (sp->xx/(1.-ConductionTask::dTemperature/mptr->pPreviousTemperature));

	// store in stress (which is minus the pressure)
	sp->xx = mPsp;
	sp->yy = mPsp;
	sp->zz = mPsp;
	
	// find the -P dV energy per unit mass dU/(rho0 V0) (uJ/g) as follows
    // dU/(rho0 V0) = - 0.5 * (pn+p(n+1))/rho0 * (V(n+1)-Vn)/V0, which simplifies to
    double dU = 0.5*(mPnsp*detf + mPsp)*(1.-1/detf);
    
    // this energy is tracked in strain energy
    mptr->AddStrainEnergy(dU);
    
    // the same energy is tracked as heat (although it will be zero if adiabatic)
    // and is dissipated (which will cause heating if adiabatic
    // Update is Cv dT - dU
    IncrementHeatEnergy(mptr,ConductionTask::dTemperature,0.,dU);
    //mptr->AddDispEnergy(dU);
        
    // the plastic energy is not otherwise used, so let's track entropy
    double dS = 0., Cv = 1000.*GetHeatCapacity(mptr);
    double Tp = mptr->pPreviousTemperature;
    double dT = ConductionTask::dTemperature;
    if(ConductionTask::energyCoupling)
    {   double dTS = dU/Cv;
        dS = Cv*(dT/Tp - dTS/(Tp+dTS));
    }
    else
        dS = (Cv*dT - dU)/Tp;
	mptr->AddPlastEnergy(dS);
}

#pragma mark IdealGas::Accessors

// return material type
const char *IdealGas::MaterialType(void) { return "Ideal Gas (Hyperelastic)"; }

// Return the material tag
int IdealGas::MaterialTag(void) { return IDEALGASMATERIAL; }

// Calculate wave speed in mm/sec.
double IdealGas::WaveSpeed(bool threeD,MPMBase *mptr)
{   return 1000.*sqrt(1.6667e9*(P0*rho)*(mptr->pTemperature/T0));
}

// calculate current wave speed in mm/sec. 
// Only change vs initial wave speed is due to J
double IdealGas::CurrentWaveSpeed(bool threeD,MPMBase *mptr)
{   // J = V/V0 = rho0/rho
    double J = mptr->GetRelativeVolume();
    return 1000.*sqrt(1.6667e9*(P0*rho/J)*(mptr->pPreviousTemperature/T0));
}





