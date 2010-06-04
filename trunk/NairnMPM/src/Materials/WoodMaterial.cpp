/********************************************************************************
	WoodMaterial.hpp
	NairnMPM

	Created by John Nairn on 6/4/10.
	Copyright 2010 Oregon State University. All rights reserved.
********************************************************************************/

#include "WoodMaterial.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark WoodMaterial::Constructors and Destructors

// Constructors
WoodMaterial::WoodMaterial() {}

/* The default contructor should call a parent class constructor and
 then fill in any new initialization.
 */
// Constructors
WoodMaterial::WoodMaterial(char *matName) : HillPlastic(matName)
{
	tempC1=100.;
	tempC2=0.;
	currentRatio=1.;
}

#pragma mark WoodMaterial::Initialization

/* If material has new property types, it must override this method and
 1. Define XML tag in the DTD file
 2. If xName matches a new property tag, set input to the type
 of variable (DOUBLE_NUM or INT_NUM) and return pointer
 to the class variable to be set.
 c. If no match, call InputMat() of superclass
 */
// Read material properties
char *WoodMaterial::InputMat(char *xName,int &input)
{
    if(strcmp(xName,"tempC1")==0)
    {	input=DOUBLE_NUM;
        return((char *)&tempC1);
    }
    else if(strcmp(xName,"tempC2")==0)
    {	input=DOUBLE_NUM;
        return((char *)&tempC2);
    }
    else
		return(HillPlastic::InputMat(xName,input));
}

/* This method is called after input file is read but before the new
 material is printed to the results file. If necessary, verify that
 all new properties are valid. If not return string with a description
 of the problem. If OK, pass on to a superclass
 NOTE: np is analysis type in case that is important
 */
// verify settings and maybe some initial calculations
const char *WoodMaterial::VerifyProperties(int np)
{
	// check at least some yielding
	if(tempC1<0.)
		return "tempC1 must be positive";
	
	// call super class
	return HillPlastic::VerifyProperties(np);
}

/* Called once at beginning (just after VerifyProperties()). For efficiency, use
 this code to calculate terms that are independent of the particle state and
 thus will remain constant throughout the calculation. Check if need to pass
 on to superclass for it to load too. Note that MaterialBase and Elastic
 set nothing so can stop if get to one of their subclasses
 */
// Constant properties used in constitutive law
void WoodMaterial::InitialLoadMechProps(int makeSpecific,int np)
{
	HillPlastic::InitialLoadMechProps(makeSpecific,np);
}

/* Print all mechanical properties or call parent class and print
 just the new mechanical properties. Use a format compatible with code
 that will read results file and similar to style of other materials
 NOTE: MaterialBase has no mechanical properties and should not be called
 NOTE: This is called after VerifyProperties() and InitialLoadMechProps()
 */
// print mechanical properties to the results
void WoodMaterial::PrintMechanicalProperties(void)
{	
    HillPlastic::PrintMechanicalProperties();
	PrintProperty("tempC1",tempC1,"");
	PrintProperty("tempC2",tempC2,"C^-1");
	cout << endl;
	
	tempC1/=100.;
	tempC2/=100.;
}

#pragma mark WoodMaterial:Methods

/* This method is called just before the constitutive law on each time step. You
 can set any parameters for that law that depends on the current state of the
 particle. Things that never change (i.e., independent of particle state) should be
 put in InitialLoadMechProps() instead.
 NOTE: For compatibility with FEA materials, some materials override LoadMechProps()
 instead, but that method only allows properties dependent on one angle (e.g., 2D
 anisotropic materials rotated about z axis).
 */
// State dependent material properties
void WoodMaterial::LoadMechanicalProps(MPMBase *mptr,int np)
{
	// scale all stiffnesses to new temperature
	double newRatio = (tempC1 + tempC2*(mptr->pPreviousTemperature-273.15));
	double newScale=newRatio/currentRatio;
	C11*=newScale;
	C12*=newScale;
	C13*=newScale;
	C22*=newScale;
	C23*=newScale;
	C33*=newScale;
	C44*=newScale;
	C55*=newScale;
	C66*=newScale;
	currentRatio=newRatio;
	
	HillPlastic::LoadMechanicalProps(mptr,np);
}

#pragma mark WoodMaterial::Custom Methods

#pragma mark WoodMaterial::Accessors

// Return the material tag
int WoodMaterial::MaterialTag(void) { return WOODMATERIAL; }

// return unique, short name for this material
const char *WoodMaterial::MaterialType(void) { return "Wood Material"; }

