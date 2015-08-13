/********************************************************************************
	WoodMaterial.cpp
	nairn-mpm-fea

	Created by John Nairn on 6/4/10.
	Copyright 2010 Oregon State University. All rights reserved.
********************************************************************************/

#include "WoodMaterial.hpp"
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
}

#pragma mark WoodMaterial::Initialization

// Read material properties
char *WoodMaterial::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    if(strcmp(xName,"tempC1")==0)
    {	input=DOUBLE_NUM;
        return((char *)&tempC1);
    }
    else if(strcmp(xName,"tempC2")==0)
    {	input=DOUBLE_NUM;
        return((char *)&tempC2);
    }
	
	return(HillPlastic::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and maybe some initial calculations
const char *WoodMaterial::VerifyAndLoadProperties(int np)
{
	if(useLargeRotation)
		return "The wood material does not yet support the large rotation options";
	
	// check at least some yielding
	if(tempC1<0.)
		return "tempC1 must be positive";
	
	// convert from percent to absolute scale
	tempC1/=100.;
	tempC2/=100.;
	
	// must call super class
	return AnisoPlasticity::VerifyAndLoadProperties(np);
}

// print mechanical properties to the results
void WoodMaterial::PrintMechanicalProperties(void) const
{	
    HillPlastic::PrintMechanicalProperties();
	PrintProperty("tempC1",tempC1*100,"");
	PrintProperty("tempC2",tempC2*100,"C^-1");
	cout << endl;
}

#pragma mark WoodMaterial:Methods

// Isotropic material can use read-only initial properties
void *WoodMaterial::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer) const
{
	AnisoPlasticProperties *p = (AnisoPlasticProperties *)AnisoPlasticity::GetCopyOfMechanicalProps(mptr,np,matBuffer,altBuffer);
	
	// calculate new ratio for reference conditions to current conditions
	double newRatio = (tempC1 + tempC2*(mptr->pPreviousTemperature-273.15));
	int i,j;
	for(i=0;i<6;i++)
		for(j=0;j<6;j++)
			p->ep.C[i][j] *= newRatio;
	
	return p;
}

#pragma mark WoodMaterial::Custom Methods

#pragma mark WoodMaterial::Accessors

// Return the material tag
int WoodMaterial::MaterialTag(void) const { return WOODMATERIAL; }

// return unique, short name for this material
const char *WoodMaterial::MaterialType(void) const { return "Wood Material"; }

