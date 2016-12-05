/******************************************************************************** 
    More MaterialBaseFEA.cpp for FEA code
    NairnFEA
    
    Created by John Nairn on Sun Mar 14 2004.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/MaterialBase.hpp"

#pragma mark MaterialBase::Initialization

// Read material properties common to all FEA materials
char *MaterialBase::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
	// but used, but allow it in case copied from MPM calculations
    if(strcmp(xName,"rho")==0)
    {	input=DOUBLE_NUM;
        return((char *)&rho);
    }

    return((char *)NULL);
}

/* calculate properties used in analyses
	If superclass overrides this method, must call this too
*/
const char *MaterialBase::VerifyAndLoadProperties(int np)
{	return NULL;
}

// print any properties common to all FEA material types
void MaterialBase::PrintCommonProperties(void) const {}

#pragma mark MaterialBase::Methods

// rotate meechanical properties if needed
void MaterialBase::LoadMechanicalPropertiesFEA(int makeSpecific,double angle,int np) { return; }

#pragma mark MaterialBase::Accessors

/* get transverse stress in plane strain analyses or get transverse strain
    in plane stress analysis
*/
double MaterialBase::GetStressStrainZZ(double sxx,double syy,double sxy,double theTemp,double angle,int np)
{	return 0.;
}

