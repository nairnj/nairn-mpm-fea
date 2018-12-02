/********************************************************************************
    RigidMaterial.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Nov 15 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.

	Rigid material in multimaterial mode
		rho = 1 mm^3/mm^3 actually, but code thinks 1 g/mm^3 or rho=1 g/mm^3
    If Planar and 3D
		mp = rho * vol (mm^3) / ptsperelement = particle volume in mm^3
        mp/rho = particle volume in mm^3 (because rho=1 g/mm^3)
    If axisymmetric
        mp = rho Ap rp / ptsperelement = particle volume per radian in mm^3
        mp/rho = particle volume per radian in mm^3 (because rho=1 g/mm^3)
        mp/(rho rp) = particle area in mm^2 (because rho=1 g/mm^3)
********************************************************************************/

#include "stdafx.h"
#include "Materials/RigidMaterial.hpp"
#include "Exceptions/CommonException.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "System/UnitsController.hpp"
#include "Read_XML/Expression.hpp"

extern double timestep;

bool RigidMaterial::someSetTemperature = false;

#pragma mark RigidMaterial::Constructors and Destructors

// Constructors with argument
RigidMaterial::RigidMaterial(char *matName,int matID,int sd) : MaterialBase(matName,matID)
{
	setDirection=sd;
	setConcentration=false;
	setTemperature=false;
	function=NULL;
	function2=NULL;
	function3=NULL;
    Vfunction=NULL;
	rho=1.;						// volume per L^3 (per mm^3 in Legacy)
	mirrored=0;
	allowsCracks=false;			// rigid material default to ignoring cracks
	useControlVelocity = false;
}

#pragma mark RigidMaterial::Initialization

// preliminary calculations (throw CommonException on problem)
const char *RigidMaterial::VerifyAndLoadProperties(int np)
{
	// is rigid multimaterial, then nothing else allowed
	if(setDirection&RIGID_MULTIMATERIAL_MODE && (setDirection!=RIGID_MULTIMATERIAL_MODE || setTemperature || setConcentration))
	{	return "Rigid material for contact in multimaterial mode cannot also set other velocities, temperature, or concentration.";
	}
	
	// force some settings for rigid contact materials
	if(setDirection&RIGID_MULTIMATERIAL_MODE)
	{	// in this mode, no boundary conditions can be applied
		setDirection = RIGID_MULTIMATERIAL_MODE;
		setConcentration = FALSE;
		setTemperature = FALSE;
		if(Vfunction!=NULL)
		{   delete Vfunction;
			Vfunction = NULL;
		}
	}
	
	// material base class does not apply to rigid materials
	//return MaterialBase::VerifyAndLoadProperties(np);
	return NULL;
}

// print to output window
void RigidMaterial::PrintMechanicalProperties(void) const
{
	char xdir='x',ydir='y';
	if(fmobj->IsAxisymmetric())
	{	xdir='R';
		ydir='Z';
	}
	
	if(setDirection&RIGID_MULTIMATERIAL_MODE)
	{	if(function!=NULL)
			cout << "Velocity " << xdir << " = " << function->GetString() << endl;
		if(function2!=NULL)
			cout << "Velocity " << ydir << " = " << function2->GetString() << endl;
		if(function3!=NULL)
			cout << "Velocity z = " << function2->GetString() << endl;
	}
	else
    {   // count velocities with functions
		int numFxns=0;
		const char *expr=NULL,*expr2=NULL,*expr3=NULL;
		if(function!=NULL)
		{   numFxns++;
			expr=function->GetString();
		}
		if(function2!=NULL)
		{   numFxns++;
			expr2=function2->GetString();
		}
		if(function3!=NULL)
		{   numFxns++;
			expr3=function3->GetString();
		}
		
        // display control directions and functinos
		if(setDirection&CONTROL_X_DIRECTION)
		{	if(setDirection&CONTROL_Y_DIRECTION)
			{	if(setDirection&CONTROL_Z_DIRECTION)
                {   // set x, y and z (with or without functions)
					cout << "Velocity in " << xdir << ", " << ydir << ", and z directions controlled" << endl;
					if(numFxns>0)
                        cout << "   Velocity " << xdir << " = " << expr << endl;
					else
						cout << "   Velocity " << xdir << " = particle velocity" << endl;
					if(numFxns>1)
                        cout << "   Velocity " << ydir << " = " << expr2 << endl;
					else
						cout << "   Velocity " << ydir << " = particle velocity" << endl;
					if(numFxns>2)
						cout << "   Velocity z = " << expr3 << endl;
					else
						cout << "   Velocity z = particle velocity" << endl;
				}
				else
                {   // set x and y (with out without function)
					cout << "Velocity in " << xdir << " and " << ydir << " directions controlled" << endl;
					if(numFxns>0)
						cout << "   Velocity " << xdir << " = " << expr << endl;
					else
						cout << "   Velocity " << xdir << " = particle velocity" << endl;
					if(numFxns>1)
						cout << "   Velocity " << ydir << " = " << expr2 << endl;
					else
						cout << "   Velocity " << ydir << " = particle velocity" << endl;
				}
			}
			else if(setDirection&CONTROL_Z_DIRECTION)
			{   // set x and z (with out without function)
				cout << "Velocity in " << xdir << " and z directions controlled" << endl;
				if(numFxns>0)
					cout << "   Velocity " << xdir << " = " << expr << endl;
				else
					cout << "   Velocity " << xdir << " = particle velocity" << endl;
				if(numFxns>1)
					cout << "   Velocity z = " << expr2 << endl;
				else
					cout << "   Velocity z = particle velocity" << endl;
			}
			else
            {   // set x direction only
				cout << "Velocity in " << xdir << " direction controlled" << endl;
                if(numFxns>0)
                    cout << "   Velocity " << xdir << " = " << expr << endl;
				else
					cout << "   Velocity " << xdir << " = particle velocity" << endl;
            }
		}
		else if(setDirection&CONTROL_Y_DIRECTION)
		{	if(setDirection&CONTROL_Z_DIRECTION)
            {   // set y and z (with or without functions)
				cout << "Velocity in " << ydir << " and z directions controlled" << endl;
				if(numFxns>0)
					cout << "   Velocity " << ydir << " = " << expr << endl;
				else
					cout << "   Velocity " << ydir << " = particle velocity" << endl;
				if(numFxns>1)
					cout << "   Velocity z = " << expr2 << endl;
				else
					cout << "   Velocity z = particle velocity" << endl;
			}
			else
            {   // set y direction only
                cout << "Velocity in " << ydir << " direction controlled" << endl;
				if(numFxns>0)
					cout << "   Velocity " << ydir << " = " << expr << endl;
				else
					cout << "   Velocity " << ydir << " = particle velocity" << endl;
            }
		}
		else if(setDirection&CONTROL_Z_DIRECTION)
        {   // set z direction only
            cout << "Velocity in z direction controlled" << endl;
			if(numFxns>0)
				cout << "   Velocity z = " << expr << endl;
			else
				cout << "   Velocity z = particle velocity" << endl;
        }

		if(!MaterialBase::extrapolateRigidBCs)
		{	if(mirrored<0)
				cout << "Velocity mirrored at minimum edge" << endl;
			else if(mirrored>0)
				cout << "Velocity mirrored at maximum edge" << endl;
		}
	}
	
	if(setTemperature || setConcentration)
    {   if(setTemperature) cout << "Temperature controlled" << endl;
        if(setConcentration)
		{	if(fmobj->HasDiffusion())
				cout << "Concentration controlled" << endl;
		}
	
        // value
		if(Vfunction!=NULL)
			cout << "Value #1 = " << Vfunction->GetString() << endl;
    }
	
	if(MaterialBase::extrapolateRigidBCs)
	{	if(setDirection&CONTROL_ANY_DIRECTION || setTemperature || setConcentration)
			cout << "Boundary conditions set after extrapolating to the grid" << endl;
	}
	
	if(setTemperature) someSetTemperature = TRUE;
	
	// optional color
	if(red>=0.)
	{	char mline[200];
		sprintf(mline,"color= %g, %g, %g, %g",red,green,blue,alpha);
		cout << mline << endl;
	}
}

// properties to read
char *RigidMaterial::InputMaterialProperty(char *xName,int &input,double &gScaling)
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
	
    else if(strcmp(xName,"SetConcentration")==0 || strcmp(xName,"SetPorePressure")==0)
    {	// also sets pore pressure if doing poroelasticity
		setConcentration=TRUE;
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
    
	else if(strcmp(xName,"mirrored")==0)
	{	input=INT_NUM;
		return((char *)&mirrored);
	}
    
    return (char *)NULL;
}

// Rigid material uses field only set to be contact material and in multimaterial mode
// throws CommonException()
int RigidMaterial::SetField(int fieldNum,bool multiMaterials,int matid,int &activeNum)
{	// not used if rigid boundary condition
	if(setDirection!=RIGID_MULTIMATERIAL_MODE)
	{	// treat as if a material that ignores cracks
		allowsCracks = true;
		return fieldNum;
	}
	
	// not allowed unless in multimaterial mode
	if(!multiMaterials)
	{	throw CommonException("Rigid material with contact not allowed in single material mode MPM",
						  "RigidMaterial::SetField");
	}
	
	// go to superclass, but do not count this rigid material as an active field
	int rigidActive=-1;			// <0 to not get pushed onto activeMatIDs
	return MaterialBase::SetField(fieldNum,multiMaterials,matid,rigidActive);
}

#pragma mark RigidMaterial::Accessors

// never called
double RigidMaterial::WaveSpeed(bool threeD,MPMBase *mptr) const { return 1.e-12; }

// return material type
const char *RigidMaterial::MaterialType(void) const
{	if(setDirection==RIGID_MULTIMATERIAL_MODE)
		return "Rigid Contact Material";
	else
		return "Rigid BC Material";
}

// return TRUE if rigid particle (for contact or for BC)
bool RigidMaterial::IsRigid(void) const { return true; }

// return TRUE if rigid BC particle (not rigid for contact)
bool RigidMaterial::IsRigidBC(void) const { return setDirection!=RIGID_MULTIMATERIAL_MODE; }

// return TRUE if rigid particle for contact
bool RigidMaterial::IsRigidContact(void) const { return setDirection==RIGID_MULTIMATERIAL_MODE; }

// return true if rigid block material
bool RigidMaterial::IsRigidBlock(void) const { return false; }

// check if should set this direction
bool RigidMaterial::RigidDirection(int aDir) const { return (setDirection&aDir)==aDir; }
bool RigidMaterial::RigidTemperature(void) const { return setTemperature; }
bool RigidMaterial::RigidConcentration(void) const { return setConcentration; }
int RigidMaterial::SetDirection(void) const { return setDirection; }

// get value function (temperature and concentration use only)
// (Don't call from parallel code due to function)
bool RigidMaterial::GetValueSetting(double *setting,double theTime,Vector *pos) const
{	if(Vfunction==NULL) return false;
	*setting=Vfunction->XYZTValue(pos,theTime*UnitsController::Scaling(1.e3));
	return true;
}

// Return true or false is this rigid particle is moving at constant velocity
// Used by ReverseLoad custom task, because only constant velocity partilces are reversed
bool RigidMaterial::IsConstantVelocity(void) const
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

// If not using control functionsm other wise evaulate functions and set
//    material point velocity in that direction (it is in vel)
// If change a direction set hasDir[i] (although really only needed for control velocity custom task)
// Note that hasDir[i] differ from set direstions only when using control velocity. If that
//     code changed to use set direction instead of controlVelocity, hasDir could be removed
// If call in parallel, use try block to catch possible function errors
bool RigidMaterial::GetVectorSetting(Vector *vel,bool *hasDir,double theTime,Vector *pos) const
{
    // false if nothing is set
    if(setDirection==0) return false;
	
	// start with false
	hasDir[0] = hasDir[1] = hasDir[2] = false;
	
    // contact rigid materials
    // false if no functions or 1 to 3 values set
	if(setDirection == RIGID_MULTIMATERIAL_MODE)
    {   if(function == NULL && function2 == NULL && function3 == NULL)
        {	// no function, but might have control velocity custom task
			if(useControlVelocity)
            {   if(controlDirection == 1)
				{	vel->x = controlVelocity;
					hasDir[0] = true;
				}
				else if (controlDirection == 2)
				{	vel->y = controlVelocity;
					hasDir[1] = true;
				}
				else if (controlDirection == 3)
				{	vel->z = controlVelocity;
					hasDir[2] = true;
				}
				else if (controlDirection == 4)
				{	*vel = ControlVelocityVector;
					hasDir[0] = true;
					hasDir[1] = true;
					hasDir[2] = true;
				}
				else
				{	// not control and no directions
					return false;
				}
				return true;
			}
			else
			{	// no function so return false
				return false;
			}
		}
		
		// here is rigid contact material and has at least one function

		// set variables
		unordered_map<string,double> vars;
		vars["t"] = theTime*UnitsController::Scaling(1.e3);
		vars["x"] = pos->x;
		vars["y"] = pos->y;
		vars["z"] = pos->z;
		vars["dt"] = timestep*UnitsController::Scaling(1.e3);
		
		// set those controlled
		if(function!=NULL)
		{	vel->x = function->EvaluateFunction(vars);
			hasDir[0] = true;
		}
		if(function2!=NULL)
		{	vel->y = function2->EvaluateFunction(vars);
			hasDir[1] = true;
		}
		if(function3!=NULL)
		{	vel->z = function2->EvaluateFunction(vars);
			hasDir[2] = true;
		}
		
		// at least on direction controlled
        return true;
	}
	
	// here will never have function(N+1) unless have function N
	
	// If no function, look for controal velocity or exit with false
	if(function == NULL)
	{	if(useControlVelocity)
		{   if(controlDirection == 1)
			{	vel->x = controlVelocity;
				hasDir[0] = true;
			}
			else if(controlDirection == 2)
			{	vel->y = controlVelocity;
				hasDir[1] = true;
			}
			else if(controlDirection == 3)
			{	vel->z = controlVelocity;
				hasDir[2] = true;
			}
			else if (controlDirection ==4)
			{	*vel = ControlVelocityVector;
				hasDir[0] = true;
				hasDir[1] = true;
				hasDir[2] = true;
			}
			else
				return false;
			return true;
		}
		return false;
	}
	
	// set variables
	unordered_map<string,double> vars;
	vars["t"] = theTime*UnitsController::Scaling(1.e3);
	vars["x"] = pos->x;
	vars["y"] = pos->y;
	vars["z"] = pos->z;
	vars["dt"] = timestep*UnitsController::Scaling(1.e3);

	// BC rigid materials - has at least one function
	// Change velocity only if have function, otherwise leave at input particle velocity
	// Note that never have functionN+1 without functionN (file reading prevents it)
	if(setDirection&CONTROL_X_DIRECTION)
	{	// x is controlled by first function
		vel->x = function->EvaluateFunction(vars);
		hasDir[0] = true;
		if(setDirection&CONTROL_Y_DIRECTION)
		{	// control x and y (change if has function #2)
			if(function2!=NULL) vel->y = function2->EvaluateFunction(vars);
			hasDir[1] = true;
			if(setDirection&CONTROL_Z_DIRECTION)
			{	// control x and y and z (change if has function #3)
				if(function3!=NULL) vel->z = function3->EvaluateFunction(vars);
				hasDir[2] = true;
			}
		}
		else if(setDirection&CONTROL_Z_DIRECTION)
		{	// control and x and z (change if has function #2)
			if(function2!=NULL) vel->z = function2->EvaluateFunction(vars);
			hasDir[2] = true;
		}
	}
	else if(setDirection&CONTROL_Y_DIRECTION)
	{	// y is controlled by first function
		vel->y = function->EvaluateFunction(vars);
		hasDir[1] = true;
		if(setDirection&CONTROL_Z_DIRECTION)
		{	// control and y and z (chagne if has function #2)
			if(function2!=NULL) vel->z = function2->EvaluateFunction(vars);
			hasDir[2] = true;
		}
	}
	else if(setDirection&CONTROL_Z_DIRECTION)
	{	// z is controlled by first function
		vel->z = function->EvaluateFunction(vars);
		hasDir[2] = true;
	}
	else
	{	// nothing was controlled (will never happen, but here for clarity)
		return false;
	}
	
	return true;
}

// setting function if needed
// throws std::bad_alloc, SAXException()
void RigidMaterial::SetSettingFunction(char *bcFunction,int functionNum)
{
	// NULL or empty is an error
	if(bcFunction==NULL)
	{	ThrowSAXException("Setting function of time and position is missing");
		return;
	}
	if(strlen(bcFunction)==0)
	{	ThrowSAXException("Setting function of time and position is missing");
		return;
	}
	
	// set one to three functions in order in check them
	switch(functionNum)
	{	case SETTING_FUNCTION_BLOCK:
			if(function!=NULL)
			{	ThrowSAXException("Duplicate setting function #1");
				return;
			}
			function =  Expression::CreateExpression(bcFunction,"Setting function #1 is not valid");
			break;
		case SETTING_FUNCTION2_BLOCK:
			if(function==NULL && setDirection!=RIGID_MULTIMATERIAL_MODE)
			{	ThrowSAXException("Cannot set function #2 before function #1 unless rigid contact material");
				return;
			}
			if(function2!=NULL)
			{	ThrowSAXException("Duplicate setting function #2");
				return;
			}
			function2 =  Expression::CreateExpression(bcFunction,"Setting function #2 is not valid");
			break;
		case SETTING_FUNCTION3_BLOCK:
			if((function==NULL || function2==NULL) && setDirection!=RIGID_MULTIMATERIAL_MODE)
			{	ThrowSAXException("Cannot set function #3 before functions #1 and #2 unless rigid contact material");
				return;
			}
			if(function3!=NULL)
			{	ThrowSAXException("Duplicate setting function #3");
				return;
			}
			function3 =  Expression::CreateExpression(bcFunction,"Setting function #3 is not valid");
			break;
        case VALUE_FUNCTION_BLOCK:
            if(Vfunction!=NULL)
			{	ThrowSAXException("Duplicate value setting function");
				return;
			}
			Vfunction =  Expression::CreateExpression(bcFunction,"Value setting function is not valid");
			break;
		default:
			break;
	}
}

// replace setting function, but HasError() is not called
void RigidMaterial::SetControlVelocity(double velocity,int direction)
{	controlVelocity = velocity;
	controlDirection = direction;
	useControlVelocity = true;
}

void RigidMaterial::SetControlVelocity(Vector *Velocity)
{	ControlVelocityVector = *Velocity;
	controlDirection = 4;
	useControlVelocity = true;
}
