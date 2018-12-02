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

// Thie method should call parent class and then fill in any new initiatizations
// Constructor
NewMaterial::NewMaterial(char *matName,int matID) : MaterialBase(matName,matID)
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

// Printing the material. The MaterialBase class prints material name and then
// Calls PrintMechanicalProperties(), PrintCommonProperties(), and PrintTransportProperties()
// MaterialBase handles common properties and transport for isotropic materials

// print mechanical properties to the results
void NewMaterial::PrintMechanicalProperties(void) const
{
	// call superclass here if it is not Material base
	
	// add new properties here
	PrintProperty("prp",newproperty,"");
	cout << endl;
}

#pragma mark NewMaterial:Step Methods

// Apply Constitutive law, check np to know what type
void NewMaterial::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res,int historyOffset) const
{
}

#pragma mark NewMaterial::Accessors

// return unique, short name for this material
const char *NewMaterial::MaterialType(void) const { return "Template Material"; }

// Calculate maximum wave speed for material in mm/sec.
double NewMaterial::WaveSpeed(bool threeD,MPMBase *mptr) const { return 1.e-12; }

