/********************************************************************************
    ImperfectInterface.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Jan 08 2006.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/ImperfectInterface.hpp"

#pragma mark ImperfectInterface::Constructors and Destructors

// Constructors
ImperfectInterface::ImperfectInterface(char *matName) : MaterialBase(matName)
{
	Dn=Dt=1.e10;
}

#pragma mark ImperfectInterface::Initialization

// Set propertties
const char *ImperfectInterface::VerifyAndLoadProperties(int np)
{
	// convert to N/m^3
	Dn*=1.e9;
	Dt*=1.e9;
	
    // Stiffness matrix
    pr.C[1][1]=Dn;
    pr.C[1][2]=Dt;
	
	return NULL;
}

// print mechanical properties to output window
void ImperfectInterface::PrintMechanicalProperties(void) const
{	
	PrintProperty("Dn",Dn/1.e9,"");
	PrintProperty("Dt",Dt/1.e9,"");
    cout << endl;
}

// Read material properties
char *ImperfectInterface::InputMat(char *xName,int &input)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"Dn")==0)
        return((char *)&Dn);
    
    else if(strcmp(xName,"Dt")==0)
        return((char *)&Dt);
    
    return MaterialBase::InputMat(xName,input);
}

#pragma mark ImperfectInterface::Accessors

// Return the material tag
int ImperfectInterface::MaterialTag(void) const { return INTERFACEPARAMS; }

// return material type
const char *ImperfectInterface::MaterialType(void) const { return "Interface parameters (in MPa/mm)"; }
