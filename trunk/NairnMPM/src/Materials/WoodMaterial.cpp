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

// Constructors
WoodMaterial::WoodMaterial(char *matName) : HillPlastic(matName)
{
	tempC1=100.;
	tempC2=0.;
	currentRatio=1.;
}

#pragma mark WoodMaterial::Initialization

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
	
	return(HillPlastic::InputMat(xName,input));
}

// verify settings and maybe some initial calculations
const char *WoodMaterial::VerifyProperties(int np)
{
	// check at least some yielding
	if(tempC1<0.)
		return "tempC1 must be positive";
	
	// must call super class
	return HillPlastic::VerifyProperties(np);
}

// Constant properties used in constitutive law
void WoodMaterial::InitialLoadMechProps(int makeSpecific,int np)
{
	HillPlastic::InitialLoadMechProps(makeSpecific,np);
}

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

// State dependent material properties
void WoodMaterial::LoadMechanicalProps(MPMBase *mptr,int np)
{
	// calculate new ratio for reference conditions to current conditions
	double newRatio = (tempC1 + tempC2*(mptr->pPreviousTemperature-273.15));
	
	// scale by new ratio (after unscaling by previos ratio)
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
	
	// remember the ratio
	currentRatio=newRatio;
	
	// call superclass
	HillPlastic::LoadMechanicalProps(mptr,np);
}

#pragma mark WoodMaterial::Custom Methods

#pragma mark WoodMaterial::Accessors

// Return the material tag
int WoodMaterial::MaterialTag(void) { return WOODMATERIAL; }

// return unique, short name for this material
const char *WoodMaterial::MaterialType(void) { return "Wood Material"; }

