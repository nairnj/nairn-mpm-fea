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

// Verify input properties do calculations; if problem return string with an error message
// If OK, MUST pass on to super class. This is called just before PrintMaterial
// (see also ValidateForUse() for checks that depend on MPM calculation mode)
const char *NewMaterial::VerifyAndLoadProperties(int np)
{
	// check properties

	// must call super class
	return MaterialBase::VerifyAndLoadProperties(np);
}

// Initialize constant transport properties used in transport calculation
//void NewMaterial::FillTransportProperties(TransportProperties *) {}

// print mechanical properties to the results
void NewMaterial::PrintMechanicalProperties(void) const
{	
	// call superclass here if it is not Material base
	
	// add new properties here
	PrintProperty("prp",newproperty,"");
    cout << endl;
}

// Print transport properties (if activated)
//void NewMaterial::PrintTransportProperties(void) const {}

// If MPM analysis not allowed, throw an exception (e.g. 3D not implemented)
// If OK, MUST pass on to super class
// (see also VerifyAndLoadProperties() for generic material property checks)
//void NewMaterial::ValidateForUse(int np) const
//{	if(np==THREED_MPM)
//	{	throw CommonException("NewMaterial cannot do 3D MPM analysis",
//							  "NewMaterial::ValidateForUse");
//	}
//	return MaterialBase::ValidateForUse(np);
//}

// If needed, a material can initialize particle state
// For example, ideal gas initializes to base line pressure
// If used, be sure to pass on to superclass when done
//void NewMaterial::SetInitialParticleState(MPMBase *mptr,int np) const {}


#pragma mark NewMaterial:HistoryVariables

// Initialize history data for a particle (if has any)
//char *NewMaterial::InitHistoryData(void) { return NULL; }

// Reutrn history data for this material type when requested (if has any)
//double NewMaterial::GetHistory(int num,char *historyPtr) const { return 0.; }

#pragma mark NewMaterial:Step Methods

// Calculate transport properties that depend on the state of the particle
//void NewMaterial::LoadTransportProps(MPMBase *mptr,int np,TransportProperties *t) const {}
//{
//}

// Implemented in case heat capacity (Cp/heat capacity for conduction) changes with particle state
// Called by conduction code
//double NewMaterial::GetHeatCapacity(MPMBase *mptr) const { return heatCapacity; }

// Apply Constitutive law, check np to know what type
void NewMaterial::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
}

#pragma mark NewMaterial::Custom Methods

#pragma mark NewMaterial::Accessors

// Return the material tag
int NewMaterial::MaterialTag(void) const { return NEWMATERIAL; }

// return unique, short name for this material
const char *NewMaterial::MaterialType(void) const { return "Template Material"; }

// Calculate maximum wave speed for material in mm/sec.
double NewMaterial::WaveSpeed(bool threeD,MPMBase *mptr) const { return 1.e-12; }

// If wave speed changes with particle state, recalculate it here and return result in mm/sec
// Only needed if wave speed changes with particle state AND if you plan to use the
// custom task to periodically adjust the time step
//double MaterialBase::CurrentWaveSpeed(bool threeD,MPMBase *mptr) const { return (new wave speed); }

// Calculate shear wave speed for material in mm/sec.
// Used only by silent boundary conditions, which are only for isotropic materials
//double NewMaterial::ShearWaveSpeed(bool threeD,MPMBase *mptr) const { return 1.e-12; }

// Maximum diffusion coefficient in cm^2/sec
//double NewMaterial::MaximumDiffusion(void) const { return 0.; }

// Maximum diffusivity in cm^2/sec
//double NewMaterial::MaximumDiffusivity(void) const { return 0.; }

// When code needs the stress, it calls this method on the material where sp is point
// to stress tensor on a particle. The default method thus returns it contents in
// a new Tensor. Some materials separated track deviatoric stress and pressure. They
// need to overide and return the full stress tensor. Other materials can use any scheme
// they want to store stress as long as this method returns the true stress
//Tensor NewMaterial::GetStress(Tensor *sp,double pressure) const { Tensor stress = *sp; return stress; }

// If new material separately track elastic and plastic strain in the particle elastic
// and plastic strain tensors, include this method and return TRUE. When true,
// code will sum the strains to get total strains, such as when finding deformation
// gradient. If FALSE, code will assume total strain is in elastic strain and therefore
// the subclass can use plastic strain for other uses (e.g. hyperelastic materials
// usual use it for the Left-Cauchy Green strain.
//bool NewMaterial::PartitionsElasticAndPlasticStrain(void) { return TRUE; }

// if a subclass material supports artificial viscosity, include this method return TRUE
// The consititutive law also needs to calculate the artificial viscosity and add
// it to pressure. A base material class method, GetArtificialViscocity(), does
// the calculation.
//bool NewMaterial::SupportsArtificialViscosity(void) const { return TRUE; }



