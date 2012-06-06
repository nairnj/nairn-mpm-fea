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
// (DOUBLE_NUM or INT_NUM) and return pointed to the class variable
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
// if OK, pass on to super class
const char *NewMaterial::VerifyProperties(int np)
{
	// check properties

	// must call super class
	return MaterialBase::VerifyProperties(np);
}

// Constant properties used in constitutive law
//void NewMaterial::InitialLoadMechProps(int makeSpecific,int np)
//{
//	MaterialBase::InitialLoadMechProps(makeSpecific,np);
//}

// Constrant transport properties used in transport calculation
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

#pragma mark NewMaterial:HistoryVariables

// Initialize history data on a particle
//char *NewMaterial::MaterialData(void) { return NULL; }

// archive history data for this material type when requested.
//double NewMaterial::GetHistory(int num,char *historyPtr) { return 0.; }

#pragma mark NewMaterial:Step Methods

// Calculate material properties that dependent on the state of the particle
//void NewMaterial::LoadMechanicalProps(MPMBase *mptr,int np)
//{
//	MaterialBase::LoadMechanicalProps(mptr,np);
//}

// Calculate transport properties that dependent on the state of the particle
//void NewMaterial::LoadTransportProps(MPMBase *mptr,int np) {}
//{
//	MaterialBase::LoadTransportProps(mptr,np);
//}

// implemented in case heat capacity (Cp or Cv) changes with particle state
//double NewMaterial::GetHeatCapacity(MPMBase *mptr) { return heatCapacity; }
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

/*
// if analysis not allowed, throw an exception
void NewMaterial::MPMConstLaw(int np)
{	if(np==THREED_MPM)
		throw CommonException("NewMaterial cannot do 3D MPM analysis","NewMaterial::MPMConstLaw");
 
	//call super class (why can't call super class?)
	return MaterialBase::MPMConstLaw(np);
}
*/

#pragma mark NewMaterial::Custom Methods

#pragma mark NewMaterial::Accessors

// Return the material tag
int NewMaterial::MaterialTag(void) { return NEWMATERIAL; }

// return unique, short name for this material
const char *NewMaterial::MaterialType(void) { return "Template Material"; }

// Calculate maximum wave speed for material in mm/sec.
double NewMaterial::WaveSpeed(bool threeD,MPMBase *mptr) { return 1.e-12; }

// Calculate shear wave speed for material in mm/sec.
//double NewMaterial::ShearWaveSpeed(bool threeD,MPMBase *mptr) { return 1.e-12; }

// maximum diffusion coefficient in cm^2/sec
//double NewMaterial::MaximumDiffusion(void) { return 0.; }

// maximum diffusivity in cm^2/sec
//double NewMaterial::MaximumDiffusivity(void) { return 0.; }


