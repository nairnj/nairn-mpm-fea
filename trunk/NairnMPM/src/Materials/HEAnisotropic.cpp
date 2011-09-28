/********************************************************************************
    HEAnisotropic.cpp
    NairnMPM
    
    Created by John Nairn, Sept 27, 2011.
    Copyright (c) 2011 John A. Nairn, All rights reserved.

	Dependencies
		HyperElastic.hpp (MaterialBase.hpp)
********************************************************************************/

#include "HEAnisotropic.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark HEAnisotropic::Constructors and Destructors

// Constructors
HEAnisotropic::HEAnisotropic() {}

/* The default contructor should call a parent class constructor and
	then fill in any new initialization.
	*/
// Constructors
HEAnisotropic::HEAnisotropic(char *matName) : HyperElastic(matName)
{
	// inintialize unique properties
}

#pragma mark HEAnisotropic::Initialization

/* This method reads properties from input file that are unique to this class.
	For details see "Creating New Material Properties" section in the
	Create_MPM_Material Wiki on nairn-fea-mem google code page
*/
// Read material properties
char *HEAnisotropic::InputMat(char *xName,int &input)
{
	// read properties for this material
    //if(strcmp(xName,"mypropertyname")==0)
    //{	input=DOUBLE_NUM;
    //    return((char *)&newproperty);
    //}
	
    return(HyperElastic::InputMat(xName,input));
}

/* This method is called after input file is read but before the new
	material is printed to the results file. If necessary, verify that
	all new properties are valid. If not return string with a description
	of the problem. If OK, must pass on to a superclass
	NOTE: np is analysis type in case that is important
*/
// verify settings and maybe some initial calculations
const char *HEAnisotropic::VerifyProperties(int np)
{
	// check properties

	// must call super class
	return HyperElastic::VerifyProperties(np);
}

/* Called once at beginning (by VerifyProperties() in MaterialBase). For efficiency,
	use this method to calculate new terms that are independent of the particle
	state and thus will remain constant throughout the calculation. When done
	(or before), pass on to super class (but MaterialBase and Elastic do not need it)
*/
// Constant properties used in constitutive law
void HEAnisotropic::InitialLoadMechProps(int makeSpecific,int np)
{
	HyperElastic::InitialLoadMechProps(makeSpecific,np);
}

/* Print all mechanical properties or call parent class and print
	just the new mechanical properties. Use a format compatible with code
	that will read results file and similar to style of other materials
	(need not pass on to MaterialBase or Elastic since they print nothing)
	NOTE: This is called after VerifyProperties() and InitialLoadMechProps()
	Sometimes scaling of properties for internal units is best done here after
	they are printed.
*/
// print mechanical properties to the results
void HEAnisotropic::PrintMechanicalProperties(void)
{	
	// call superclass here if it is not Material base
	HyperElastic::PrintMechanicalProperties();
	
	// add new properties here
}

#pragma mark HEAnisotropic:Methods

/*	Apply 2D constitutive law updating all needed terms for material type. Required updates are:
		stress, strain, plastic strain (all components) (stress should be a specific stress)
		rotation strain (single angle)
		strain energy, plastic energy, and dissipated energy (dissipated needed if want to couple to conduction)
	To support thermal and solvent expansion, include their effect on strains
    If there are material-related data on the particle, update them too
	dvij are (gradient rates X time increment) to give deformation gradient change
*/
void HEAnisotropic::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
        double delTime,int np)
{
}

/* Apply 3D constitutive law updating all needed terms for material type. Required updates are:
		stress, strain, plastic strain (all components) (stress should be a specific stress)
		rotation strain
		strain energy, plastic energy, and dissipated energy (dissipated needed if want to couple to conduction)
	To support thermal and solvent expansion, include their effect on strains
    If there are material-related data on the particle, update them too
    dvij are (gradient rates X time increment) to give deformation gradient change
*/
void HEAnisotropic::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
        double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np)
{
}

#pragma mark HEAnisotropic::Custom Methods

#pragma mark HEAnisotropic::Accessors

// Return the material tag
int HEAnisotropic::MaterialTag(void) { return HEANISOTROPIC; }

// return unique, short name for this material
const char *HEAnisotropic::MaterialType(void) { return "Hyperelastic Anisotropic"; }

// If this material supports 3D MPM, then remove this method
bool HEAnisotropic::ThreeDMaterial(void) { return true; }

/* Calculate maximum wave speed for material in mm/sec. WaveSpeed called only
	once for each material point at beginning of calculation. If variable wave
	speed, be conservative and return the maximum possible save speed.
*/
// wave speed for this material
double HEAnisotropic::WaveSpeed(bool threeD) { return 1.e-12; }

/* This method is only used by silent boundary conditions. The MaterialBase base
	class returns WaveSpeed()/sqrt(3). Override only if have better result
*/
// shear wave speed for this material
//double HEAnisotropic::ShearWaveSpeed(bool threeD) { return 1.e-12; }

