/********************************************************************************
    BoundaryCondition.cpp
    NairnMPM
    
    Created by John Nairn on Mon Mar 29 2004.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include "Read_XML/mathexpr.hpp"

// global expression variables
double BoundaryCondition::varTime=0.;
double BoundaryCondition::varXValue=0.;
double BoundaryCondition::varYValue=0.;
double BoundaryCondition::varZValue=0.;
double BoundaryCondition::varRotValue=0.;
PRVar varTimeArray[5] = { NULL, NULL, NULL, NULL, NULL };

#pragma mark BoundaryCondition: Constructors and Destructors

BoundaryCondition::BoundaryCondition(int bcStyle,double bcValue,double bcTime)
{
    style=bcStyle;
    value=bcValue;
    ftime=bcTime;
	offset=0.;
	function=NULL;
	scale=1.;				// for function settings
}

// Reuse Rigid properties (subclass set other properties) and return this for calling function use
BoundaryCondition *BoundaryCondition::SetRigidProperties(long num,int dof,int bcStyle,double bcValue)
{	nodeNum=num;
    style=bcStyle;
    value=bcValue;
	ftime=0.;
	return (BoundaryCondition *)this;
}

// Destructor (and it is virtual)
BoundaryCondition::~BoundaryCondition()
{	if(function!=NULL) delete function;
}

// just unset condition, because may want to reuse it, return next one to unset
BoundaryCondition *BoundaryCondition::UnsetDirection(void)
{	// return next one
	return (BoundaryCondition *)GetNextObject();
}

#pragma mark BoundaryCondition: Methods

// Calculate value for time in milleseconds
double BoundaryCondition::BCValue(double mstime)
{
	double currentValue=0.;
	
    switch(style)
    {   case CONSTANT_VALUE:
            if(mstime>=ftime) currentValue=value;
            break;
        case LINEAR_VALUE:
            if(mstime>=ftime) currentValue=offset+value*(mstime-ftime);
            break;
        case SINE_VALUE:
			currentValue=value*sin(mstime*ftime);
            break;
        case COSINE_VALUE:
            currentValue=value*cos(mstime*ftime);
            break;
		case FUNCTION_VALUE:
			if(mstime>=ftime)
			{	varTime=mstime-ftime;
				GetPosition(&varXValue,&varYValue,&varZValue,&varRotValue);
				currentValue=scale*function->Val();
			}
        default:
            break;
    }
    return currentValue;
}

// print function if needed and a new line
void BoundaryCondition::PrintFunction(ostream &os)
{	if(style==FUNCTION_VALUE)
	{	char *expr=GetFunctionString();
		os << "  " << expr;
		delete [] expr;
	}
	os << endl;
}

#pragma mark BoundaryCondition: Accessors

// if BC is activated, return nodeNum, otherwise return 0
long BoundaryCondition::GetNodeNum(double mstime)
{
    switch(style)
    {   case CONSTANT_VALUE:
        case LINEAR_VALUE:
            if(mstime<ftime) return 0L;
            break;
        case SINE_VALUE:
        case COSINE_VALUE:
        default:
            break;
    }
	return nodeNum;
}

// just return the nodeNum
long BoundaryCondition::GetNodeNum(void) { return nodeNum; }

// set function if needed
void BoundaryCondition::SetFunction(char *bcFunction)
{
	if(style!=FUNCTION_VALUE) return;
	if(bcFunction==NULL)
		throw SAXException("Boundary condition function of time and position is missing");
	if(strlen(bcFunction)==0)
		throw SAXException("Boundary condition function of time and position is missing");
	if(function!=NULL)
		throw SAXException("Duplicate boundary condition functions of time");
	
	// create variable
	if(varTimeArray[0]==NULL)
	{	varTimeArray[0]=new RVar("t",&varTime);
		varTimeArray[1]=new RVar("x",&varXValue);
		varTimeArray[2]=new RVar("y",&varYValue);
		varTimeArray[3]=new RVar("z",&varZValue);
		varTimeArray[4]=new RVar("q",&varRotValue);
	}
		
	// create the function and check it
	function=new ROperation(bcFunction,5,varTimeArray);
	if(function->HasError())
		throw SAXException("Boundary condition function of time and position is not valid");
	
	// keep value, which can option by + or - to tell visualization code the direction of the load
	// value=1.;
}

// new string for the function (caller should delete it)
char *BoundaryCondition::GetFunctionString(void)
{	if(function!=NULL) return function->Expr('#');
	char *unknown=new char[2];
	unknown[0]='?';
	unknown[1]=0;
	return unknown;
}

// get position for current boundary conditions
// assumes nodal BC, override for particle BC
void BoundaryCondition::GetPosition(double *xpos,double *ypos,double *zpos,double *rot)
{	int i=GetNodeNum();
	*xpos=nd[i]->x;
	*ypos=nd[i]->y;
	*zpos=nd[i]->z;
	*rot=0.;
}


