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
	
	// Ideal gas has heat capacity in J/(kg-K)
	// Must set C to Cv to work well with theory, even though output and equations
	//    for conductivity claim Cp is needed. For Ideal Gas, Cv = 1.5R
	// This material overrides any attempt to set heat capacity
	heatCapacity = 1500.*P0/(T0*rho);
	heatCapacityVol = heatCapacity;
	
    // call super class
    return MaterialBase::VerifyProperties(np);
}

// if analysis not allowed, throw an exception
void IdealGas::ValidateForUse(int np)
{	if(np==PLANE_STRESS_MPM)
	{	throw CommonException("IdealGas material cannot do 2D plane stree MPM analysis",
							  "NairnMPM::ValidateForUse");
	}
	
	// call super class (why can't call super class?)
	return MaterialBase::ValidateForUse(np);
}

// Private properties used in constitutive law
void IdealGas::InitialLoadMechProps(int makeSpecific,int np)
{
	hasMatProps=TRUE;
	
	// P0 in specific units for MPM of N/m^2 cm^3/g
	P0sp=P0*1.0e+06/rho;
	
}

// If needed, a material can initialize particle state
// For example, ideal gas initializes to base line pressure
void IdealGas::SetInitialParticleState(MPMBase *mptr,int np)
{
    // Find Cauchy stresses
    double Psp = P0sp * (thermal.reference/T0);
	Tensor *sp=mptr->GetStressTensor();
	sp->xx = -Psp;
	sp->yy = -Psp;
	sp->zz = -Psp;
}


#pragma mark IdealGas::Methods

/* For 2D MPM analysis, take increments in strain and calculate new
    Particle: strains, rotation strain, stresses, strain energy, angle
    dvij are (gradient rates X time increment) to give deformation gradient change
	Does not support thermal or moisture strains
*/
void IdealGas::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
								double delTime,int np)
{
	// get new deformation gradient
	double F[3][3];
	double detf = GetDeformationGrad(F,mptr,dvxx,dvyy,dvxy,dvyx,TRUE,TRUE);

	// compute specific pressure p/rho0 = P0sp * (T/T0) * (1/J)
	// incrementally J(n+1) = det f Jn, p(n+1)/rho0 = P0sp * (T(n+1)/T0) * (1/(det f Jn))
	//												= P0sp * (Tn/T0) * (1/Jn) * (T(n+1)/Tn) * (1/det f)
	//												= pn * (T(n+1)/Tn) * (1/det f)
	//double Psp = - P0sp * (mptr->pPreviousTemperature/T0) / J;
	Tensor *sp=mptr->GetStressTensor();
	double Psp = (sp->xx/(1.-ConductionTask::dTemperature/mptr->pPreviousTemperature)) / detf;
	
	// find (Cauchy stress)/rho0 (if CONSTANT_RHO) or (Kirchoff stress)/rho0 (if not)
	sp->xx = Psp;
	sp->yy = Psp;
	sp->zz = Psp;
	
	// internal energy change (per rho0) = - 0.5 * (pn+p(n+1))/rho0 * (V(n+1)-Vn)
	double dU = -0.5*(P0sp/T0)*(mptr->pPreviousTemperature*(detf-1./detf)
								- ConductionTask::dTemperature*(detf-1.));
	mptr->AddStrainEnergy(dU);
	mptr->AddDispEnergy(dU);
}

/* For 3D MPM analysis, take increments in strain and calculate new
	Particle: strains, rotation strain, stresses, strain energy, angle
	dvij are (gradient rates X time increment) to give deformation gradient change
*/
void IdealGas::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
						  double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np)
{
	// get determinent of incremental deformation gradient (and update strains)
	double F[3][3];
	double detf = GetDeformationGrad(F,mptr,dvxx,dvyy,dvzz,dvxy,dvyx,dvxz,dvzx,dvyz,dvzy,TRUE,TRUE);

	// compute specific pressure p/rho0 = P0sp * (T/T0) * (1/J)
	// incrementally J(n+1) = det f Jn, p(n+1)/rho0 = P0sp * (T(n+1)/T0) * (1/(det f Jn))
	//												= P0sp * (Tn/T0) * (1/Jn) * (T(n+1)/Tn) * (1/det f)
	//												= pn * (T(n+1)/Tn) * (1/det f)
	//double Psp = - P0sp * (mptr->pPreviousTemperature/T0) / J;
	Tensor *sp=mptr->GetStressTensor();
	double Psp = (sp->xx/(1.-ConductionTask::dTemperature/mptr->pPreviousTemperature)) / detf;
	
	// Find Cauchy stresses
	sp->xx = Psp;
	sp->yy = Psp;
	sp->zz = Psp;
   
	// internal energy change (per rho0) = - 0.5 * (pn+p(n+1))/rho0 * (V(n+1)-Vn)
	double dU = -0.5*(P0sp/T0)*(mptr->pPreviousTemperature*(detf-1./detf)
										- ConductionTask::dTemperature*(detf-1.));
	mptr->AddStrainEnergy(dU);
	mptr->AddDispEnergy(dU);
}

#pragma mark IdealGas::Accessors

// return material type
const char *IdealGas::MaterialType(void) { return "Ideal Gas (Hyperelastic)"; }

// Return the material tag
int IdealGas::MaterialTag(void) { return IDEALGASMATERIAL; }

// Calculate wave speed in mm/sec. Here using adiabatic bulk modulus at the starting temperature
// If T rises a lot during compression, the problem may become unstable. The solution is
// to anticipate and use some time step saftelt factor. The problem will probably not arise
// when gas using in conjuctions with solids having much higher wave speeds.
double IdealGas::WaveSpeed(bool threeD,MPMBase *mptr)
{	return sqrt(1.6667e9*(P0/rho)*(mptr->pTemperature/T0));
}




