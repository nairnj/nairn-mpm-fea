/********************************************************************************
 IdealGas.cpp
 NairnMPM
 
 Created by Edward Le and Peter Mackenzie on May 26, 2012.
 Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/IdealGas.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
 
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
    if(P0 <= 0. || rho <= 0.0 || T0 <= 0.0 )
		return "Ideal gas material model needs positive parameters P0, rho, and T0";

    // set current temperature to initial temperature 
    // until a better value becomes available 
    // (stored as a member variable)
    Temp = T0;

    // call super class
    return MaterialBase::VerifyProperties(np);
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
	// current temperature in Kelvin (now stored as a member variable)
	Temp = mptr->pPreviousTemperature;
	
	// get new deformation gradient
	double F[3][3];
	double detf = GetDeformationGrad(F,mptr,dvxx,dvyy,dvxy,dvyx,TRUE,TRUE);

	// left Cauchy deformation tensor B = F F^T
	Tensor B = GetLeftCauchyTensor2D(F);
	
	//double J2 = B.xx*B.yy*B.zz + 2.*B.xy*B.xz*B.yz - B.yz*B.yz*B.xx - B.xz*B.xz*B.yy - B.xy*B.xy*B.zz;
	//double J = sqrt(J2);
	/*
	John,

	The above step is not needed since detf == J (if everything is computed correctly).

	Peter
	*/
	
	// get new Jacobian determinant and update strains and rotations
	double J = 1.0;
#ifndef CONSTANT_RHO
	// J = GetCurrentRelativeVolume(mptr,FALSE);	
        J = detf;
#endif
	
	/*
	 John,
	 
	 A thought on adiabatic simulations: 
	   o temperature can be computed from the adiabatic equation of state
	     instead of a thermal analysis and plugged in the next line.  This would do adiabatic changes, 
	     but ignores any thermal calculation in the background.
	   o instead, a heat source could be added to the thermal analysis and its magniture would need
	     to be computed from within the material class.  This would allow for non-equilibrium computations that 
	     include limits of isothermal and adiabatic changes of state.  Don't know if that is of interest
	     to you, or even if that's possible in the current code.
	*/
	
	// compute specific pressure p/rho0
	//double Psp = 0.0;	
	double Psp = P0sp;	// plane stress default: set to ambient pressure at T0
	if(np==PLANE_STRAIN_MPM)
	{	
		Psp = P0sp * (Temp/T0) / J;
	}
	
	// find (Cauchy stress)/rho0 (if CONSTANT_RHO) or (Kirchoff stress)/rho0 (if not)
	Tensor *sp=mptr->GetStressTensor();
	sp->xx = -Psp;
	sp->yy = -Psp;
	sp->xy =  0.0;
	if(np==PLANE_STRAIN_MPM)
	{	sp->zz =  -Psp;
	}
	
	// strain energy (total energy divided by initial rho)
	double energy;
	energy = J*(J*J*J-1.0)*Psp/3.0;
	mptr->SetStrainEnergy(energy);
}

/* For 3D MPM analysis, take increments in strain and calculate new
	Particle: strains, rotation strain, stresses, strain energy, angle
	dvij are (gradient rates X time increment) to give deformation gradient change
*/
void IdealGas::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
						  double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np)
{
	// current temperature in Kelvin (now stored as a member variable)
	Temp = mptr->pPreviousTemperature;

	// get new deformation gradient
	double F[3][3];
	double detf = GetDeformationGrad(F,mptr,dvxx,dvyy,dvzz,dvxy,dvyx,dvxz,dvzx,dvyz,dvzy,TRUE,TRUE);

        double J = 1.0;
#ifndef CONSTANT_RHO
        double J = detf;
#endif

	// left Cauchy deformation tensor B = F F^T
	Tensor B = GetLeftCauchyTensor3D(F);
	
	//double J2 = B.xx*B.yy*B.zz + 2.*B.xy*B.xz*B.yz - B.yz*B.yz*B.xx - B.xz*B.xz*B.yy - B.xy*B.xy*B.zz;
	//double J = sqrt(J2);
	/*
	John,

	The above step is not needed since detf == J (if everything is computed correctly).

	Peter
	*/
	
	// compute specific pressure p/rho0
	double Psp = P0sp * (Temp/T0) / J;
	
	// Find Cauchy stresses
	Tensor *sp=mptr->GetStressTensor();
	sp->xx = -Psp;
	sp->yy = -Psp;
	sp->zz = -Psp;
	sp->xy = 0.0;
	sp->xz = 0.0;
	sp->yz = 0.0;
    
	// strain energy
	double energy;
	energy = J*(J*J*J-1.0)*Psp/3.0;
	mptr->SetStrainEnergy(energy);
}

#pragma mark IdealGas::Accessors

// return material type
const char *IdealGas::MaterialType(void) { return "Ideal Gas (Hyperelastic)"; }

// Return the material tag
int IdealGas::MaterialTag(void) { return IDEALGASMATERIAL; }

//	calculate wave speed in mm/sec
double IdealGas::WaveSpeed(bool threeD)
{
    return sqrt(1.e9*(P0/rho)*(Temp/T0));
}




