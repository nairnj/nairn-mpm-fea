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
double RigidMaterial::xPos=0.;
double RigidMaterial::yPos=0.;
double RigidMaterial::zPos=0.;
PRVar rmTimeArray[4] = { NULL,NULL,NULL,NULL };

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
	function2=NULL;
	function3=NULL;
    Vfunction=NULL;
	rho=1000.;
}

#pragma mark RigidMaterial::Initialization

// preliminary calculations (throw CommonException on problem)
const char *RigidMaterial::VerifyProperties(int np)
{	// is rigid multimaterial, then nothing else allowed
	if(setDirection&RIGID_MULTIMATERIAL_MODE && (setDirection!=RIGID_MULTIMATERIAL_MODE || setTemperature || setConcentration))
	{	return "Rigid material for contact in multimaterial mode cannot also set other velocities, temperature, or concentration.";
	}
	
	// to super class
	return MaterialBase::VerifyProperties(np);
}

// print to output window
void RigidMaterial::PrintMechanicalProperties(void)
{
	if(setDirection&RIGID_MULTIMATERIAL_MODE)
	{	cout << "Rigid multimaterial with contact" << endl;
        if(function!=NULL)
        {	char *expr=function->Expr('#');
			cout << "Velocity x = " << expr << endl;
            delete [] expr;
        }
        if(function2!=NULL)
        {	char *expr=function2->Expr('#');
			cout << "Velocity y = " << expr << endl;
            delete [] expr;
        }
        if(function3!=NULL)
        {	char *expr=function3->Expr('#');
			cout << "Velocity z = " << expr << endl;
            delete [] expr;
        }
        
		// in this mode, no boundary conditions can be applied
		setDirection=RIGID_MULTIMATERIAL_MODE;
		setConcentration=FALSE;
		setTemperature=FALSE;
        if(Vfunction!=NULL)
        {   delete Vfunction;
            Vfunction=NULL;
        }
	}
	else
    {   // count velocities with functions
		int numFxns=0;
        char *expr=NULL,*expr2=NULL,*expr3=NULL;
		if(function!=NULL)
        {   numFxns++;
            expr=function->Expr('#');
        }
		if(function2!=NULL)
        {   numFxns++;
            expr2=function2->Expr('#');
        }
		if(function3!=NULL)
        {   numFxns++;
            expr3=function3->Expr('#');
        }
        
        // display control directions and functinos
		if(setDirection&CONTROL_X_DIRECTION)
		{	if(setDirection&CONTROL_Y_DIRECTION)
			{	if(setDirection&CONTROL_Z_DIRECTION)
                {   // set x, y and z (unless too few functions)
					if(numFxns==1)
                    {   cout << "Velocity in x direction controlled" << endl;
                        cout << "Velocity x = " << expr << endl;
                    }
					else if(numFxns==2)
                    {   cout << "Velocity in x and y directions controlled" << endl;
                        cout << "Velocity x = " << expr << endl;
                        cout << "Velocity y = " << expr2 << endl;
                    }
					else
                    {   cout << "Velocity in x, y, and z directions controlled" << endl;
                        if(numFxns>2)
                        {   cout << "Velocity x = " << expr << endl;
                            cout << "Velocity y = " << expr2 << endl;
                            cout << "Velocity z = " << expr3 << endl;
                        }
                    }
				}
				else
                {   // set x and y (unless to few functions)
					if(numFxns==1)
                    {   cout << "Velocity in x direction controlled" << endl;
                        cout << "Velocity x = " << expr << endl;
                    }
					else
                    {   cout << "Velocity in x and y directions controlled" << endl;
                        if(numFxns>1)
                        {   cout << "Velocity x = " << expr << endl;
                            cout << "Velocity y = " << expr2 << endl;
                        }
                    }
				}
			}
			else if(setDirection&CONTROL_Z_DIRECTION)
            {   // set x and z (unless too few functions)
				if(numFxns==1)
                {   cout << "Velocity in x direction controlled" << endl;
                    cout << "Velocity x = " << expr << endl;
                }
				else
                {   cout << "Velocity in x and z directions controlled" << endl;
                    if(numFxns>1)
                    {   cout << "Velocity x = " << expr << endl;
                        cout << "Velocity z = " << expr2 << endl;
                    }
                }
			}
			else
            {   // set x direction only
				cout << "Velocity in x direction controlled" << endl;
                if(numFxns>=1)
                    cout << "Velocity x = " << expr << endl;
            }
		}
		else if(setDirection&CONTROL_Y_DIRECTION)
		{	if(setDirection&CONTROL_Z_DIRECTION)
            {   // set y and z (unless too few functions)
				if(numFxns==1)
                {   cout << "Velocity in y direction controlled" << endl;
                    cout << "Velocity y = " << expr << endl;
                }
				else
                {   cout << "Velocity in y and z directions controlled" << endl;
                    if(numFxns>1)
                    {   cout << "Velocity y = " << expr << endl;
                        cout << "Velocity z = " << expr2 << endl;
                    }
                }
			}
			else
            {   // set y direction only
                cout << "Velocity in y direction controlled" << endl;
                if(numFxns>=1)
                    cout << "Velocity y = " << expr << endl;
            }
		}
		else if(setDirection&CONTROL_Z_DIRECTION)
        {   // set z direction only
            cout << "Velocity in z direction controlled" << endl;
            if(numFxns>=1)
                cout << "Velocity z = " << expr << endl;
        }
        
        // clean up
        if(expr!=NULL) delete [] expr;
        if(expr2!=NULL) delete [] expr2;
        if(expr3!=NULL) delete [] expr3;
	}
	
	if(setTemperature || setConcentration)
    {   if(setTemperature) cout << "Temperature controlled" << endl;
        if(setConcentration) cout << "Concentration controlled" << endl;
	
        // value
        if(Vfunction!=NULL)
        {	char *expr=Vfunction->Expr('#');
			cout << "Value #1 = " << expr << endl;
            delete [] expr;
        }
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
	
	else if(strcmp(xName,"SettingFunction2")==0)
	{	input=SETTING_FUNCTION2_BLOCK;
		return((char *)this);
	}
	
	else if(strcmp(xName,"SettingFunction3")==0)
	{	input=SETTING_FUNCTION3_BLOCK;
		return((char *)this);
	}
    
	else if(strcmp(xName,"ValueFunction")==0)
	{	input=VALUE_FUNCTION_BLOCK;
		return((char *)this);
	}
    
    return (char *)NULL;
}

// Rigid material uses field only set to be contact material and in multimaterial mode
int RigidMaterial::SetField(int fieldNum,bool multiMaterials,int matid)
{	// not used if rigid bpundary condition
	if(setDirection!=RIGID_MULTIMATERIAL_MODE) return fieldNum;
	
	// not allowed unless in multimaterial mode
	if(!multiMaterials)
	{	throw CommonException("Rigid material with contact not allowed in single material mode MPM",
						  "RigidMaterial::SetField");
	}
	
	// go to superclass
	return MaterialBase::SetField(fieldNum,multiMaterials,matid);
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
double RigidMaterial::WaveSpeed(bool threeD,MPMBase *mptr) { return 1.e-12; }

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

// get value function (temperature and concentration use only)
bool RigidMaterial::GetValueSetting(double *setting,double theTime,Vector *pos)
{	if(Vfunction==NULL) return false;
	varTime=1000.*theTime;
	xPos=pos->x;
	yPos=pos->y;
	zPos=pos->z;
	*setting=Vfunction->Val();
	return true;
}

// Return true or false is this rigid particle is moving at constant velocity
// Used by ReverseLoad custom task, because only constant velocity partilces are reversed
bool RigidMaterial::IsConstantVelocity(void)
{
    // If no velocity is set, perhaps temperature or concentration are being set
    if(setDirection==0)
	{	if(Vfunction==NULL) return true;
		return false;
	}
	
    // contact rigid materials: true if no functions or 1 to 3 values set
	// rigid BCs: true is no first function
	if(setDirection==RIGID_MULTIMATERIAL_MODE)
    {   if(function==NULL && function2==NULL && function3==NULL)
			return true;
	}
	else if(function==NULL)
		return true;
	
	// the particle motion is controlled by a function
	return false;
}

// get vector from one to three functions and the directions being set in hasDir
// If no directions set or no functions in place return false
bool RigidMaterial::GetVectorSetting(Vector *vel,bool *hasDir,double theTime,Vector *pos)
{
    // false if nothing is set
    if(setDirection==0) return false;
	
    // contact rigid materials
    // false if no functions or 1 to 3 values set
	if(setDirection==RIGID_MULTIMATERIAL_MODE)
    {   if(function==NULL && function2==NULL && function3==NULL)
            return false;
        
        // set variables
        varTime=1000.*theTime;
        xPos=pos->x;
        yPos=pos->y;
        zPos=pos->z;
	
		if(function!=NULL)
        {	vel->x=function->Val();
            hasDir[0]=true;
        }
		if(function2!=NULL)
		{	vel->y = function2->Val();
			hasDir[1]=true;
		}
		if(function3!=NULL)
		{	vel->z = function3->Val();
			hasDir[2]=true;
		}
        
        return true;
	}
    
    // if no functions, then done
    if(function==NULL) return false;
    
	// set one to three components (in order of functions and set directions)
	hasDir[0] = hasDir[1] = hasDir[2] = false;
    
    // set variables
    varTime=1000.*theTime;
    xPos=pos->x;
    yPos=pos->y;
    zPos=pos->z;
    
    // BC rigid materials
	if(setDirection&CONTROL_X_DIRECTION)
	{	vel->x=function->Val();
		hasDir[0]=true;
		if(function2!=NULL)
		{	if(setDirection&CONTROL_Y_DIRECTION)
			{	vel->y = function2->Val();
				hasDir[1]=true;
				if(function3!=NULL)
				{	if(setDirection&CONTROL_Z_DIRECTION)
					{	vel->z = function3->Val();
						hasDir[2]=true;
					}
				}
			}
			else if(setDirection&CONTROL_Z_DIRECTION)
			{	vel->z = function2->Val();
				hasDir[2]=true;
			}
		}
	}
	else if(setDirection&CONTROL_Y_DIRECTION)
	{	vel->y=function->Val();
		hasDir[1]=true;
		if(function2!=NULL)
		{	if(setDirection&CONTROL_Z_DIRECTION)
			{	vel->z = function2->Val();
				hasDir[2]=true;
			}
		}
	}
	else if(setDirection&CONTROL_Z_DIRECTION)
	{	vel->z=function->Val();
		hasDir[2]=true;
	}
    
	return true;
}

// setting function if needed
void RigidMaterial::SetSettingFunction(char *bcFunction,int functionNum)
{
	// NULL or empty is an error
	if(bcFunction==NULL)
		ThrowSAXException("Setting function of time and position is missing");
	if(strlen(bcFunction)==0)
		ThrowSAXException("Setting function of time and position is missing");
	
	// create time variable if needed
	if(rmTimeArray[0]==NULL)
	{	rmTimeArray[0]=new RVar("t",&varTime);
		rmTimeArray[1]=new RVar("x",&xPos);
		rmTimeArray[2]=new RVar("y",&yPos);
		rmTimeArray[3]=new RVar("z",&zPos);
	}
	
	// set one to three functions in order in check them
	switch(functionNum)
	{	case SETTING_FUNCTION_BLOCK:
			if(function!=NULL)
				ThrowSAXException("Duplicate setting function #1");
			function=new ROperation(bcFunction,4,rmTimeArray);
			if(function->HasError())
				ThrowSAXException("Setting function #1is not valid");
			break;
		case SETTING_FUNCTION2_BLOCK:
			if(function==NULL && setDirection!=RIGID_MULTIMATERIAL_MODE)
				ThrowSAXException("Cannot set function #2 before function #1 unless rigid contact material");
			if(function2!=NULL)
				ThrowSAXException("Duplicate setting function #2");
			function2=new ROperation(bcFunction,4,rmTimeArray);
			if(function2->HasError())
				ThrowSAXException("Setting function #2 is not valid");
			break;
		case SETTING_FUNCTION3_BLOCK:
			if((function==NULL || function2==NULL) && setDirection!=RIGID_MULTIMATERIAL_MODE)
				ThrowSAXException("Cannot set function #3 before functions #1 and #2 unless rigid contact material");
			if(function3!=NULL)
				ThrowSAXException("Duplicate setting function #3");
			function3=new ROperation(bcFunction,4,rmTimeArray);
			if(function3->HasError())
				ThrowSAXException("Setting function #3 is not valid");
			break;
        case VALUE_FUNCTION_BLOCK:
            if(Vfunction!=NULL)
				ThrowSAXException("Duplicate value setting function");
			Vfunction=new ROperation(bcFunction,4,rmTimeArray);
			if(Vfunction->HasError())
				ThrowSAXException("Value setting function is not valid");
			break;
		default:
			break;
	}
	
}


