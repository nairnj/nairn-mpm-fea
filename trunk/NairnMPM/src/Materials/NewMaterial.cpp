/********************************************************************************
    NewMaterial.cpp
    NairnMPM
    
    Created by John Nairn, June 16, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		Orthotropic.hpp (TranIsotropic.hpp, MaterialBase.hpp)
********************************************************************************/

#include "NewMaterial.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark NewMaterial::Constructors and Destructors

// Constructors
NewMaterial::NewMaterial() {}

/* The default contructor should call a parent class constructor and
	then fill in any new initialization.
	*/
// Constructors
NewMaterial::NewMaterial(char *matName) : MaterialBase(matName)
{
	newproperty=0.;
}

#pragma mark NewMaterial::Initialization

/* If material has new property types, it must override this method and
	1. Define XML tag in the DTD file
	2. If xName matches a new property tag, set input to the type
		of variable (DOUBLE_NUM or INT_NUM) and return pointer
		to the class variable to be set.
	c. If no match, call InputMat() of superclass
*/
// Read material properties
char *NewMaterial::InputMat(char *xName,int &input)
{
    if(strcmp(xName,"newproperty")==0)
    {	input=DOUBLE_NUM;
        return((char *)&newproperty);
    }
    
    return(MaterialBase::InputMat(xName,input));
}

/* This method is called after input file is read but before the new
	material is printed to the results file. If necessary, verify that
	all new properties are valid. If not return string with a description
	of the problem. If OK, pass on to a superclass
	NOTE: np is analysis type in case that is important
*/
// verify settings and maybe some initial calculations
const char *NewMaterial::VerifyProperties(int np)
{
	// check properties

	// call super class
	return MaterialBase::VerifyProperties(np);
}

/* This method called before starting, but only if this material is actually
	assigned to at least one material point. If it cannot be used in the
	current analysis type, throw CommonException("error message","code"). If OK,
	call ValidateUse() in a superclass
*/
// if cannot be used in current analysis type throw MPMTermination()
void NewMaterial::ValidateUse(int np)
{
	MaterialBase::ValidateUse(np);
}

/* Called once at beginning (just after VerifyProperties()). For efficiency, use
	this code to calculate terms that are independent of the particle state and
	thus will remain constant throughout the calculation. Check if need to pass
	on to superclass for it to load too. Note that MaterialBase and Elastic
	set nothing so can stop if get to one of their subclasses
*/
// Constant properties used in constitutive law
void NewMaterial::InitialLoadMechProps(int makeSpecific,int np) { }

/* Print all mechanical properties or call parent class and print
	just the new mechanical properties. Use a format compatible with code
	that will read results file and similar to style of other materials
	NOTE: MaterialBase has no mechanical properties and should not be called
	NOTE: This is called after VerifyProperties() and InitialLoadMechProps()
*/
// print mechanical properties to the results
void NewMaterial::PrintMechanicalProperties(void)
{	
	PrintProperty("prp",newproperty,"");
    cout << endl;
}

/* This method only needed if this material needs history dependent properties. If so,
	create space for them (it is called once for each particle with that material),
	initialize, and return pointer. If override similar method in parent class,
	will need to make sure all history data is correctly allocated and accessible.
*/
// Initialize history data
//char *NewMaterial::MaterialData(void) { return NULL; }

/* Print all transport properties in format compatible with code that will read
	results file. Only print if transprot is activated
	NOTE: Base class prints isotropic properties if transport is activated. If this
	is enough, no additional printing is needed.
*/
// Print transport properties
//void NewMaterial::PrintTransportProperties(void) {}

/* Called once at beginning (just after VerifyProperties()). For efficiency, use
	this code to calculate transport terms that are independent of the particle state
	and thus will remain constant throughout the calculation
	NOTE: The base material class automatically supports isotropic diffusion and conduction
	using material properties D (in diffusionCon) for diffusion constant and
	kCond (in kCond) for thermal conductivity. If this is enough, no additional
	work is needed.
	*/
// Constrant transport properties used in transport calculation
//void NewMaterial::InitialLoadTransProps(void) {}

#pragma mark NewMaterial:Methods

/* This method is called just before the constitutive law on each time step. You
	can set any parameters for that law that depends on the current state of the
	particle. Things that never change (i.e., independent of particle state) should be
	put in InitialLoadMechProps() instead.
	NOTE: For compatibility with FEA materials, some materials override LoadMechProps()
		instead, but that method only allows properties dependent on one angle (e.g., 2D
		anisotropic materials rotated about z axis).
*/
// State dependent material properties
//void NewMaterial::LoadMechanicalProps(MPMBase *mptr,int np) {}

/* Called when looping over material points to store parameters needed in transport calculations (tensors).
	It is called prior to transport task to AddForces(). Only needed for anistropic materials
	or those whose transport properties change depending on particle state. Load the properties into diffusionTensor
	and kCondTensor in the MaterialBase class. Things that never change (i.e., independent of particle state)
	should be put in InitialLoadTransProps() instead.
*/
// State dependent material properties
//void NewMaterial::LoadTransportProps(MPMBase *mptr) {}

/* When conduction is activated, this method is called before calculations that depend
	on heat capacity. If it changes with particle state, return new result in units
	of J/(g-K).
	NOTE: The heatCapacity and heatCapacityVol properties in base class are converted to
		these units when the analysis starts.
	NOTE: GetHeatCapacityVol() for Cv not used by any materials, but could be implemented
		and called if needed.
*/
// implemented in case heat capacity (Cp or Cv) changes with particle state
//double MaterialBase::GetHeatCapacity(MPMBase *mptr) { return heatCapacity; }
//double MaterialBase::GetHeatCapacityVol(MPMBase *mptr) { return heatCapacityVol; }

/*	Apply 2D constitutive law updating all needed terms for material type Required updates are:
		stress, strain, plastic strain (all components) (stress should be a specific stress)
		rotation strain (single angle)
		strain energy, plastic energy, and dissipated energy (dissipated needed if want to couple to conduction)
	To support thermal and solvent expansion, include their effect on strains
    If there are material-related data on the particle, update them too
	dvij are (gradient rates X time increment) to give deformation gradient change
*/
void NewMaterial::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
        double delTime,int np)
{
}

/* Apply 3D constitutive law updating all needed terms for material type Required updates are:
		stress, strain, plastic strain (all components) (stress should be a specific stress)
		rotation strain
		strain energy, plastic energy, and dissipated energy (dissipated needed if want to couple to conduction)
	To support thermal and solvent expansion, include their effect on strains
    If there are material-related data on the particle, update them too
    dvij are (gradient rates X time increment) to give deformation gradient change
*/
void NewMaterial::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
        double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np)
{
}

#pragma mark NewMaterial::Cutom Methods

#pragma mark NewMaterial::Accessors

// Return the material tag
int NewMaterial::MaterialTag(void) { return NEWMATERIAL; }

// return unique, short name for this material
const char *NewMaterial::MaterialType(void) { return "Template Material"; }

// If this material supports 3D MPM, then remove this method
bool NewMaterial::ThreeDMaterial(void) { return false; }

/* If the material puts history variables on the particle, this method should
    return history variable number num stored in the material data pointed to by
    historyPtr. If num is invalid return 0. It is only used to archive history
	data and the archiving can currently only archive history #1
*/
// archive history data for this material type when requested.
//double NewMaterial::GetHistory(int num,char *historyPtr) { return 0.; }

/* Calculate maximum wave speed for material in mm/sec. WaveSpeed called only
	once for each material point at beginning of calculation. If variable wave
	speed, be conservative and return the maximum possible save speed.
*/
// wave speed for this material
double NewMaterial::WaveSpeed(bool threeD) { return 1.e-12; }

/* This method is only used by silent boundary conditions. The MaterialBase base
	class returns WaveSpeed()/sqrt(3). Override only if have better result
*/
// shear wave speed for this material
//double NewMaterial::ShearWaveSpeed(bool threeD) { return 1.e-12; }

/* Calculate maximum diffusion constant in cm^2/sec. Method called once for each material
	point at beginning of calculation
	NOTE: MaterialBase handles isotropic diffusion. No need to override if that is
	enough
*/
// maximum diffusion coefficient in cm^2/sec
//double NewMaterial::MaximumDiffusion(void) { return 0.; }

/* Calculate maximum diffusivity in cm^2/sec (= k/(100 rho Cp)). Method called once for 
	each material point at beginning of calculation
	NOTE: MaterialBase handles isotropic conduction. No need to override if that is
	enough
*/
// maximum diffusivity in cm^2/sec
//double NewMaterial::MaximumDiffusivity(void) { return 0.; }

