/********************************************************************************
    NewMaterial.cpp
    NairnMPM
    
    Created by John Nairn, June 16, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
 
	See Create_MPM_Material for details in google code web site wiki

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#include "NewMaterial.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark NewMaterial::Constructors and Destructors

// Constructors
NewMaterial::NewMaterial() {}

// The default contructor should call a parent class constructor and
// then fill in any new initialization.
NewMaterial::NewMaterial(char *matName) : MaterialBase(matName)
{
	newproperty=0.;
}

#pragma mark NewMaterial::Initialization

// Read material properties by name (in xName). Set input to variable type
// (DOUBLE_NUM or INT_NUM) and return pointer to the class variable
// (cast as a char *)
char *NewMaterial::InputMat(char *xName,int &input)
{
	// read properties for this material
    if(strcmp(xName,"newproperty")==0)
    {	input=DOUBLE_NUM;
        return((char *)&newproperty);
    }
	
    return(MaterialBase::InputMat(xName,input));
}

// Verify input properties; if problem return string with an error message
// If OK, MUST pass on to super class
// (see also ValidateForUse() for checkes that depend on MPM calculation mode)
const char *NewMaterial::VerifyProperties(int np)
{
	// check properties

	// must call super class
	return MaterialBase::VerifyProperties(np);
}

// Initialize constant properties used in constitutive law
//void NewMaterial::InitialLoadMechProps(int makeSpecific,int np)
//{
//	MaterialBase::InitialLoadMechProps(makeSpecific,np);
//}

// Initialize constant transport properties used in transport calculation
//void NewMaterial::InitialLoadTransProps(void) {}

// print mechanical properties to the results
void NewMaterial::PrintMechanicalProperties(void)
{	
	// call superclass here if it is not Material base
	
	// add new properties here
	PrintProperty("prp",newproperty,"");
    cout << endl;
}

// Print transport properties (if activated)
//void NewMaterial::PrintTransportProperties(void) {}

// If MPM analysis not allowed, throw an exception (e.g. 3D not implemented)
// If OK, MUST pass on to super class
// (see also VerifyProperties() for generic material property checks)
//void NewMaterial::ValidateForUse(int np)
//{	if(np==THREED_MPM)
//	{	throw CommonException("NewMaterial cannot do 3D MPM analysis",
//							  "NewMaterial::ValidateForUse");
//	}
//	return MaterialBase::ValidateForUse(np);
//}

#pragma mark NewMaterial:HistoryVariables

// Initialize history data on a particle (if has any)
//char *NewMaterial::MaterialData(void) { return NULL; }

// Reutrn history data for this material type when requested (if has any)
//double NewMaterial::GetHistory(int num,char *historyPtr) { return 0.; }

#pragma mark NewMaterial:Step Methods

// Calculate material properties that depend on the state of the particle
// If implemented, MUST pass onto super class
//void NewMaterial::LoadMechanicalProps(MPMBase *mptr,int np)
//{
//	MaterialBase::LoadMechanicalProps(mptr,np);
//}

// Calculate transport properties that depend on the state of the particle
// If implemented, MUST pass onto super class
//void NewMaterial::LoadTransportProps(MPMBase *mptr,int np) {}
//{
//	MaterialBase::LoadTransportProps(mptr,np);
//}

// Implemented in case heat capacity (Cp/heat capacity for conduction) changes with particle state
// Called by conduction code
//double NewMaterial::GetHeatCapacity(MPMBase *mptr) { return heatCapacity; }

// Implemented in case heat capacity Cv changes with particle state (but not called
//		by any core methods)
//double NewMaterial::GetHeatCapacityVol(MPMBase *mptr) { return heatCapacityVol; }

// Apply 2D Constitutive law
void NewMaterial::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
        double delTime,int np)
{
}

// Apply 3D Constitutive law
void NewMaterial::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
        double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np)
{
}

#pragma mark NewMaterial::Custom Methods

#pragma mark NewMaterial::Accessors

// Return the material tag
int NewMaterial::MaterialTag(void) { return NEWMATERIAL; }

// return unique, short name for this material
const char *NewMaterial::MaterialType(void) { return "Template Material"; }

// Calculate maximum wave speed for material in mm/sec.
double NewMaterial::WaveSpeed(bool threeD,MPMBase *mptr) { return 1.e-12; }

// If wave speed changes with particle state, recalculate it here and return result in mm/sec
// Only needed if wave speed changes with particle state AND if you plan to use the
// custom task to periodically adjust the time step
//double MaterialBase::CurrentWaveSpeed(bool threeD,MPMBase *mptr) { return (new wave speed); }

// Calculate shear wave speed for material in mm/sec.
// Used only by silent boundary conditions, which are only for isotropic materials
//double NewMaterial::ShearWaveSpeed(bool threeD,MPMBase *mptr) { return 1.e-12; }

// Maximum diffusion coefficient in cm^2/sec
//double NewMaterial::MaximumDiffusion(void) { return 0.; }

// Maximum diffusivity in cm^2/sec
//double NewMaterial::MaximumDiffusivity(void) { return 0.; }


