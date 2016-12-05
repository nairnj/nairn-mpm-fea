/********************************************************************************
    ImperfectInterface.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Jan 08 2006.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/ImperfectInterface.hpp"
#include "System/UnitsController.hpp"

#pragma mark ImperfectInterface::Constructors and Destructors

// Constructors
ImperfectInterface::ImperfectInterface(char *matName) : MaterialBase(matName)
{
	Dn=Dt=1.e15;
}

#pragma mark ImperfectInterface::Initialization

// Set propertties
const char *ImperfectInterface::VerifyAndLoadProperties(int np)
{
    // Stiffness matrix
    pr.C[1][1]=Dn;
    pr.C[1][2]=Dt;
	
	return NULL;
}

// print mechanical properties to output window
void ImperfectInterface::PrintMechanicalProperties(void) const
{	
	PrintProperty("Dn",Dn,UnitsController::Label(INTERFACEPARAM_UNITS));
	PrintProperty("Dt",Dt,UnitsController::Label(INTERFACEPARAM_UNITS));
    cout << endl;
}

// Read material properties
char *ImperfectInterface::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"Dn")==0)
        return((char *)&Dn);
    
    else if(strcmp(xName,"Dt")==0)
        return((char *)&Dt);
    
    return MaterialBase::InputMaterialProperty(xName,input,gScaling);
}

#pragma mark ImperfectInterface::Accessors

// return material type
const char *ImperfectInterface::MaterialType(void) const { return "Interface parameters"; }
