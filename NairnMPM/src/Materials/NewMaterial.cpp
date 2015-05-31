/********************************************************************************
    NewMaterial.cpp
    nairn-mpm-fea
    
    Created by John Nairn, June 16, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
 
	See Create_MPM_Material for details in google code web site wiki
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
char *NewMaterial::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
	// read properties for this material
    if(strcmp(xName,"newproperty")==0)
    {	input=DOUBLE_NUM;
        return((char *)&newproperty);
    }
	
    return(MaterialBase::InputMaterialProperty(xName,input,gScaling));
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

// if analysis not allowed, throw a CommonException
// Call super class when checks are done
//bool NewMaterial::ValidateForUse(int np) const
//{	if(np==THREED_MPM)
//	{	throw CommonException("NewMaterial cannot do 3D MPM analysis",
//							  "NewMaterial::ValidateForUse");
//	}
//	MaterialBase::ValidateForUse(np);
//}

#pragma mark NewMaterial:HistoryVariables

// If needed, a material can initialize particle state
// For example, ideal gas initializes to base line pressure
// If used, be sure to pass on to superclass when done
//void NewMaterial::SetInitialParticleState(MPMBase *mptr,int np) const {}

// Initialize history data for a particle (if has any)
//char *NewMaterial::InitHistoryData(void) { return NULL; }

// Reutrn history data for this material type when requested (if has any)
//double NewMaterial::GetHistory(int num,char *historyPtr) const { return 0.; }

#pragma mark NewMaterial:Step Methods

// Return buffer size and buffer size for hardenling law (if used)
//int MaterialBase::SizeOfMechanicalProperties(int &altBufferSize) const
//{   return 0;
//}

// Get copy of properties in the material class that depend on material state
//void *MaterialBase::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer) const
//{	return NULL;
//}

// Calculate transport properties that depend on the state of the particle
//void NewMaterial::LoadTransportProps(MPMBase *mptr,int np,TransportProperties *t) const {}
//{
//}

// Implemented in case heat capacity (Cp/heat capacity for conduction) changes with particle state
// Called by conduction code and return nJ/(g-K)
//double NewMaterial::GetHeatCapacity(MPMBase *mptr) const { return heatCapacity; }

// A material can override to set Cp-Cv in nJ/(g-K)
// mptr will be NULL when printing material properties
//double NewMaterial::GetCpMinusCv(MPMBase *mptr) const { return 0; }

// Apply Constitutive law, check np to know what type
void NewMaterial::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
}

#pragma mark NewMaterial::Custom Methods

#pragma mark NewMaterial::Accessors

// return unique, short name for this material
const char *NewMaterial::MaterialType(void) const { return "Template Material"; }

// Return the material tag
int NewMaterial::MaterialTag(void) const { return NEWMATERIAL; }

// Calculate maximum wave speed for material in mm/sec.
double NewMaterial::WaveSpeed(bool threeD,MPMBase *mptr) const { return 1.e-12; }

// If wave speed changes with particle state, recalculate it here and return result in mm/sec
// Only needed if wave speed changes with particle state AND if you plan to use the
// custom task to periodically adjust the time step
//double MaterialBase::CurrentWaveSpeed(bool threeD,MPMBase *mptr) const { return (new wave speed); }

// Calculate shear wave speed for material in mm/sec.
// Used only by silent boundary conditions, which are only for isotropic materials
//double NewMaterial::ShearWaveSpeed(bool threeD,MPMBase *mptr) const { return 1.e-12; }

// Maximum diffusion coefficient in mm^2/sec
//double NewMaterial::MaximumDiffusion(void) const { return 0.; }

// Maximum diffusivity in mm^2/sec
//double NewMaterial::MaximumDiffusivity(void) const { return 0.; }

// When code needs the stress, it calls this method on the material where sp is point
// to stress tensor on a particle. The default method thus returns it contents in
// a new Tensor. Some materials separated track deviatoric stress and pressure. They
// need to overide and return the full stress tensor. Other materials can use any scheme
// they want to store stress as long as this method returns the true stress
//Tensor NewMaterial::GetStress(Tensor *sp,double pressure) const { Tensor stress = *sp; return stress; }

// if a subclass material supports artificial viscosity, include this method return TRUE
// The consititutive law also needs to calculate the artificial viscosity and add
// it to pressure. A base material class method, GetArtificialViscocity(), does
// the calculation.
//bool NewMaterial::SupportsArtificialViscosity(void) const { return TRUE; }



