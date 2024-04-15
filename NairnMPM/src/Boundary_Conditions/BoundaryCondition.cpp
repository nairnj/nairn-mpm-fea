/********************************************************************************
    BoundaryCondition.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Mon Mar 29 2004.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/
#ifdef _MSWIN_
#include "stdafx.h"
#endif 
#include <cmath>

#include "Boundary_Conditions/BoundaryCondition.hpp"

#define MPM_CODE
#include "Nodes/NodalPoint.hpp"
#include "System/UnitsController.hpp"
#include "Exceptions/CommonException.hpp"
#include "Read_XML/Expression.hpp"

#pragma mark BoundaryCondition: Constructors and Destructors

BoundaryCondition::BoundaryCondition(int bcStyle,double bcValue,double bcTime)
{
    style = bcStyle;
    SetBCValue(bcValue);
	SetBCFirstTime(bcTime);
	
	offset = 0.;
	function = NULL;
	scale = 1.;				// scale function calculation only if needed
	bcID = 0;
}

// Reuse Rigid properties (subclass set other properties) and return this for calling function use
// Set using internal units so no need for legacy scaling and always using CONSTANT_VALUE
BoundaryCondition *BoundaryCondition::SetRigidProperties(int num,int dof,int bcStyle,double bcValue)
{	// set direction
    nd[num]->SetFixedDirection(dof);
    
    // set other properties
    nodeNum = num;
    style = bcStyle;
    value = bcValue;
	ftime = 0.;
	bcID = 0;
	return (BoundaryCondition *)this;
}

// Destructor (and it is virtual)
BoundaryCondition::~BoundaryCondition()
{	if(function!=NULL) delete function;
}

// get set direction
int BoundaryCondition::GetSetDirection(void) const { return 0; }

// just unset condition, because may want to reuse it, return next one to unset
BoundaryCondition *BoundaryCondition::UnsetDirection(void)
{	nd[nodeNum]->UnsetFixedDirection(GetSetDirection());
	return (BoundaryCondition *)GetNextObject();
}

#pragma mark BoundaryCondition: Methods

// print it (if can be used)
BoundaryCondition *BoundaryCondition::PrintBC(ostream &os)
{
    char nline[200];
    snprintf(nline,200,"%7d %2d %15.7e %15.7e",nodeNum,style,GetBCValueOut(),GetBCFirstTimeOut());
    os << nline;
    PrintFunction(os);
    return (BoundaryCondition *)GetNextObject();
}

// Calculate value at stepTime
double BoundaryCondition::BCValue(double stepTime)
{
	double currentValue=0.;
	
    switch(style)
    {   case CONSTANT_VALUE:
            if(stepTime>=ftime) currentValue=value;
            break;
        case LINEAR_VALUE:
            if(stepTime>=ftime) currentValue=offset+value*(stepTime-ftime);
            break;
        case SINE_VALUE:
			currentValue=value*sin(stepTime*ftime);
            break;
        case COSINE_VALUE:
            currentValue=value*cos(stepTime*ftime);
            break;
		case FUNCTION_VALUE:
			if(stepTime>=ftime)
			{
                // (see Expression vmap)
				double vars[7];
				vars[0] = 6.5;
				vars[1] = UnitsController::Scaling(1.e3)*(stepTime-ftime);		//t
				GetPositionVars(vars);
				currentValue = scale*function->EvaluateFunction(vars);
			}
        default:
            break;
    }
 
    return currentValue;
}

// print function if needed and a new line
void BoundaryCondition::PrintFunction(ostream &os)
{	if(style==FUNCTION_VALUE)
	{	if(function!=NULL)
			os << " " << function->GetString();
		else
			os << " " << "?";
	}
	os << endl;
}

#pragma mark BoundaryCondition: Accessors

// if BC is activated, return nodeNum, otherwise return 0
int BoundaryCondition::GetNodeNum(double bctime)
{
    switch(style)
    {   case CONSTANT_VALUE:
        case LINEAR_VALUE:
            if(bctime<ftime) return 0L;
            break;
        case SINE_VALUE:
        case COSINE_VALUE:
        default:
            break;
    }
    return nd[nodeNum]->NodeHasParticles() ? nodeNum : 0L ;
}

// just return the nodeNum
int BoundaryCondition::GetNodeNum(void) { return nodeNum; }

// set function if needed
// throws std::bad_alloc, SAXException()
void BoundaryCondition::SetFunction(char *bcFunction)
{
	// skip if not needed
	if(style!=FUNCTION_VALUE) return;
	
	// if needed, cannot be NULL, empty, or a duplicate
	if(bcFunction==NULL)
	{	ThrowSAXException("Boundary condition function of time and position is missing");
		return;
	}
	if(strlen(bcFunction)==0)
	{	ThrowSAXException("Boundary condition function of time and position is missing");
		return;
	}
	if(function!=NULL)
	{	ThrowSAXException("Duplicate boundary condition functions of time");
		return;
	}

	function =  Expression::CreateExpression(bcFunction,"Boundary condition function not valid");
	
	// keep value, which can option by + or - to tell visualization code the direction of the load
	// value=1.;
}

// put x,y,z,q into fixed places in array (see Expression vmap)
void BoundaryCondition::GetPositionVars(double *vars)
{	int i=GetNodeNum();
	vars[2] = nd[i]->x;		//x
	vars[3] = nd[i]->y;		//y
	vars[4] = nd[i]->z;		//z
	vars[6] = 0.;			//q
}

// Boundary condition ID may be used for some purpose by certain conditions
// Be sure to set it whenever create or reuse a boundary conditions
int BoundaryCondition::GetID(void) { return bcID; }
void BoundaryCondition::SetID(int newID) { bcID = newID; }

// Accessors

double *BoundaryCondition::GetBCValuePtr(void) { return &value; }
void BoundaryCondition::SetBCValue(double bcvalue)
{	// Legacy linear value in value/ms need to be converted to value/sec
	if(style==LINEAR_VALUE)
		value = bcvalue*UnitsController::Scaling(1.e3);
	else
		value = bcvalue;
}
void BoundaryCondition::SetBCValueCU(double bcvalue) { value = bcvalue; }
double BoundaryCondition::GetBCValue(void) { return value; }
double BoundaryCondition::GetBCValueOut(void)
{	// Legacy linear value in value/s need to be converted to value/ms for legacy output only
	if(style==LINEAR_VALUE)
		return UnitsController::Scaling(1.e-3)*value;
	else
		return value;
}

// check for 1/time in sine and cosine styles
void BoundaryCondition::SetBCFirstTime(double bcftime)
{	// Legacy frequencies 1/ms to 1/s and times ms to s
	if(style==SINE_VALUE || style==COSINE_VALUE)
		ftime = UnitsController::Scaling(1.e3)*bcftime;
	else
		ftime = UnitsController::Scaling(1.e-3)*bcftime;
}
void BoundaryCondition::SetBCFirstTimeCU(double bcftime) { ftime = bcftime; }
double BoundaryCondition::GetBCFirstTime(void) { return ftime; }
double BoundaryCondition::GetBCFirstTimeOut(void)
{	// Legacy frequencies 1/s to 1/ms and times s to ms for legacy output only
	if(style==SINE_VALUE || style==COSINE_VALUE)
		return UnitsController::Scaling(1.e-3)*ftime;
	else
		return UnitsController::Scaling(1.e3)*ftime;
}

void BoundaryCondition::SetBCOffset(double bcoffset) { offset = bcoffset; }
double BoundaryCondition::GetBCOffset(void) { return offset; }
int BoundaryCondition::GetBCStyle(void) { return style; }
