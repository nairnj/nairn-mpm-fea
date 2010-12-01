/********************************************************************************
    RigidMaterial.cpp
    NairnMPM
    
    Created by John Nairn on Nov 15 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.

	Rigid material in multimaterial mode
		rho = 1000 mm^3/cm^3
		mp = rho * vol (cm^3) / ptsperelement = particle volume in mm^3
********************************************************************************/

#include "Materials/RigidMaterial.hpp"
#include "Read_XML/mathexpr.hpp"
#include "Exceptions/CommonException.hpp"

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
	rho=1000.;
}

#pragma mark RigidMaterial::Initialization

// print to output window
void RigidMaterial::PrintMechanicalProperties(void)
{
	if(setDirection&RIGID_MULTIMATERIAL_MODE)
	{	cout << "Rigid multimaterial with contact" << endl;
		
		// in this mode, no boundary conditions can be applied
		setDirection=RIGID_MULTIMATERIAL_MODE;
		setConcentration=FALSE;
		setTemperature=FALSE;
	}
	else
	{	if(setDirection&CONTROL_X_DIRECTION && setDirection&CONTROL_Y_DIRECTION)
			cout << "Velocity in skewed x-y direction controlled" << endl;
		else if(setDirection&CONTROL_X_DIRECTION)
			cout << "Velocity in x direction controlled" << endl;
		else if(setDirection&CONTROL_Y_DIRECTION)
			cout << "Velocity in y direction controlled" << endl;
		if(setDirection&CONTROL_Z_DIRECTION)
			cout << "Velocity in z direction controlled" << endl;
	}
	if(setTemperature) cout << "Temperature controlled" << endl;
	if(setConcentration) cout << "Concentration controlled" << endl;
	if(function!=NULL)
	{	char *expr=function->Expr('#');
		if(setDirection&RIGID_MULTIMATERIAL_MODE)
			cout << "Velocity = ";
		else
			cout << "Value = ";
		cout << expr << endl;
		delete [] expr;
	}
	
	// optional color
	if(red>=0.)
	{	char mline[200];
		sprintf(mline,"color= %g, %g, %g",red,green,blue);
		cout << mline << endl;
	}
}

// properties to read
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

// Rigid material uses field only set to be contact material and in multimaterial mode
int RigidMaterial::SetField(int fieldNum,bool multiMaterials,int matid)
{	// not used if rigid bpundary condition
	if(setDirection!=RIGID_MULTIMATERIAL_MODE) return fieldNum;
	
	// not allowed unless in multimaterial mode
	if(!multiMaterials) return -1;
	
	// go to superclass
	return MaterialBase::SetField(fieldNum,multiMaterials,matid);
}

// preliminary calculations (throw CommonException on problem)
void RigidMaterial::PreliminaryMatCalcs(void)
{	// is rigid multimaterial, then nothing else allowed
	if(setDirection&RIGID_MULTIMATERIAL_MODE && (setDirection!=RIGID_MULTIMATERIAL_MODE || setTemperature || setConcentration))
		throw CommonException("Rigid material for contact in multimaterial mode cannot also set other velocities, temperature, or concentration.","RigidMaterial::PreliminaryMatCalcs");
	MaterialBase::PreliminaryMatCalcs();
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
double RigidMaterial::WaveSpeed(bool threeD) { return 1.e-12; }

// Return the material tag
int RigidMaterial::MaterialTag(void) { return RIGIDMATERIAL; }

// return material type
const char *RigidMaterial::MaterialType(void) { return "Rigid Material"; }

// return TRUE if rigid particle (for contact or for BC)
short RigidMaterial::Rigid(void) { return TRUE; }

// return TRUE if rigid BC particle (not rigid for contact)
short RigidMaterial::RigidBC(void) { return setDirection!=RIGID_MULTIMATERIAL_MODE; }

// return TRUE if rigid particle for contact
short RigidMaterial::RigidContact(void) { return setDirection==RIGID_MULTIMATERIAL_MODE; }

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
		ThrowSAXException("Setting function of time is missing");
	if(strlen(bcFunction)==0)
		ThrowSAXException("Setting function of time is missing");
	if(function!=NULL)
		ThrowSAXException("Duplicate setting functions of time");
	
	// create variable
	if(rmTimeArray[0]==NULL)
	{	rmTimeArray[0]=new RVar("t",&varTime);
	}
		
	// create the function and check it
	function=new ROperation(bcFunction,1,rmTimeArray);
	if(function->HasError())
		ThrowSAXException("Setting function of time is not valid");
}


