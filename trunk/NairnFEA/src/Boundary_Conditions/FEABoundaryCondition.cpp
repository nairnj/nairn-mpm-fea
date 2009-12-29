/********************************************************************************
    FEABoundaryCondition.cpp
    NairnFEA
    
    Created by John Nairn on July 24, 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		mathexpr
********************************************************************************/

#include "Read_XML/CommonReadHandler.hpp"
#include "FEABoundaryCondition.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Read_XML/mathexpr.hpp"

// global expression variables
double FEABoundaryCondition::varXValue=0.;
double FEABoundaryCondition::varYValue=0.;
double FEABoundaryCondition::varRValue=0.;
double FEABoundaryCondition::varZValue=0.;
PRVar varArray[4] = { NULL, NULL, NULL, NULL };

#pragma mark FEABoundaryCondition: Constructors and Destructors

FEABoundaryCondition::FEABoundaryCondition()
{
}

// Destructor (and it is virtual)
FEABoundaryCondition::~FEABoundaryCondition()
{	
}

#pragma mark FEABoundaryCondition: Methods

// decrement node number or delete this BC if no longer user
// deleted will be true only if first object was deleved
LinkedObject *FEABoundaryCondition::DecrementNodeNum(long removedNode,LinkedObject *prevObject,bool *deleted)
{
	*deleted=FALSE;
	
	// if uses this node, remove from the list and delete it
	if(nodeNum==removedNode)
	{	*deleted=TRUE;
		if(prevObject!=NULL)
			prevObject->SetNextObject(nextObject);
		LinkedObject *returnObject=nextObject;
		delete this;
		return returnObject;
	}
	
	// just decrement if higher node number
	else if(nodeNum>removedNode)
		nodeNum--;
	
	return nextObject;
}
		
#pragma mark FEABoundaryCondition: Accessors

void FEABoundaryCondition::SetValue(double numValue,char *bcFunction)
{
	// take number if no value
	if(bcFunction==NULL)
	{	bcValue=numValue;
		return;
	}
	
	// decode function at nodal position
	if(strlen(bcFunction)==0)
		throw SAXException("Boundary condition function of position is empty");
	
	// create variable
	if(varArray[0]==NULL)
	{	varArray[0]=new RVar("x",&varXValue);
		varArray[1]=new RVar("y",&varYValue);
		varArray[2]=new RVar("r",&varRValue);
		varArray[3]=new RVar("z",&varZValue);
	}
		
	// create the function and check it
	ROperation *function=new ROperation(bcFunction,4,varArray);
	if(function->HasError())
		throw SAXException("Boundary condition function of position is not valid");
	
	// set value to 1, it might be used to scale the function result
	varXValue=varRValue=nd[nodeNum]->x;
	varYValue=varZValue=nd[nodeNum]->y;
	bcValue=function->Val();
	delete function;
}

// return pointer to the bcValue - a double
double *FEABoundaryCondition::GetValuePtr(void) { return &bcValue; }

