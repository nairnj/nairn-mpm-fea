/********************************************************************************
    RigidMaterial.cpp
    NairnMPM
    
    Created by John Nairn on Nov 15 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.

	Dependencies
********************************************************************************/

#include "Materials/RigidMaterial.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include "Read_XML/mathexpr.hpp"

#pragma mark RigidMaterial::Constructors and Destructors

// global expression variables
double RigidMaterial::varTime=0.;
PRVar rmTimeArray[4] = { NULL };

// Constructors
RigidMaterial::RigidMaterial()
{
}

// Constructors with arguments 
RigidMaterial::RigidMaterial(char *matName) : MaterialBase(matName)
{
	setDirection=0;
	setConcentration=FALSE;
	setTemperature=FALSE;
	function=NULL;
}

#pragma mark RigidMaterial::Initialization

// print to output window
void RigidMaterial::PrintMechanicalProperties(void)
{
	if(setDirection&1 && setDirection&2)
		cout << "Velocity in skewed x-y direction controlled" << endl;
	else if(setDirection&1)
		cout << "Velocity in x direction controlled" << endl;
	else if(setDirection&2)
		cout << "Velocity in y direction controlled" << endl;
	if(setDirection&4)
		cout << "Velocity in z direction controlled" << endl;
	if(setTemperature) cout << "Temperature controlled" << endl;
	if(setConcentration) cout << "Concentration controlled" << endl;
	if(function!=NULL)
	{	char *expr=function->Expr('#');
		cout << "Value = " << expr << endl;
		delete [] expr;
	}
}

// no properties to read
char *RigidMaterial::InputMat(char *xName,int &input)
{
    if(strcmp(xName,"SetDirection")==0)
    {	input=INT_NUM;
        return((char *)&setDirection);
    }
	
    else if(strcmp(xName,"SetTemperature")==0)
    {	setTemperature=TRUE;
        input=NOT_NUM;
        return((char *)&setTemperature);
    }
	
    else if(strcmp(xName,"SetConcentration")==0)
    {	setConcentration=TRUE;
        input=NOT_NUM;
        return((char *)&setConcentration);
    }

	else if(strcmp(xName,"SettingFunction")==0)
	{	input=SETTING_FUNCTION_BLOCK;
		return((char *)this);
	}
	
    return (char *)NULL;
}

#pragma mark RigidMaterial::Methods

// does not have 2D constutuve law
void RigidMaterial::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvxy,double dvyx,
        double delTime,int np)
{}

// Does not have 3D constitutive law
void RigidMaterial::MPMConstLaw(MPMBase *mptr,double dvxx,double dvyy,double dvzz,double dvxy,double dvyx,
        double dvxz,double dvzx,double dvyz,double dvzy,double delTime,int np)
{}

#pragma mark RigidMaterial::Accessors

// never called
double RigidMaterial::WaveSpeed(void) { return 1.e-12; }

// Return the material tag
int RigidMaterial::MaterialTag(void) { return RIGIDMATERIAL; }

// return material type
const char *RigidMaterial::MaterialType(void) { return "Rigid Material"; }

// check if rigid material
short RigidMaterial::Rigid(void) { return TRUE; }

// check if should set this direction
bool RigidMaterial::RigidDirection(int aDir) { return (setDirection&aDir)==aDir; }
bool RigidMaterial::RigidTemperature(void) { return setTemperature; }
bool RigidMaterial::RigidConcentration(void) { return setConcentration; }

// get scaling factor if a function is defined
bool RigidMaterial::GetSetting(double *setting,double theTime)
{	if(function==NULL) return false;
	varTime=1000.*theTime;
	*setting=function->Val();
	return true;
}

// set function if needed
void RigidMaterial::SetSettingFunction(char *bcFunction)
{
	if(bcFunction==NULL)
		throw SAXException("Setting function of time is missing");
	if(strlen(bcFunction)==0)
		throw SAXException("Setting function of time is missing");
	if(function!=NULL)
		throw SAXException("Duplicate Setting functions of time");
	
	// create variable
	if(rmTimeArray[0]==NULL)
	{	rmTimeArray[0]=new RVar("t",&varTime);
	}
		
	// create the function and check it
	function=new ROperation(bcFunction,1,rmTimeArray);
	if(function->HasError())
		throw SAXException("Setting function of time is not valid");
}


