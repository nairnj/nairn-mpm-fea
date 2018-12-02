/********************************************************************************
    FEABoundaryCondition.cpp
    NairnFEA
    
    Created by John Nairn on July 24, 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		mathexpr
********************************************************************************/

#include "stdafx.h"
#include "FEABoundaryCondition.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Read_XML/Expression.hpp"

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
LinkedObject *FEABoundaryCondition::DecrementNodeNum(int removedNode,LinkedObject *prevObject,bool *deleted)
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

// set the value
// throws std::bad_alloc, SAXException()
void FEABoundaryCondition::SetValue(double numValue,char *bcFunction)
{
	// take number if no value
	if(bcFunction==NULL)
	{	bcValue=numValue;
		return;
	}
	
	// decode function at nodal position
	if(strlen(bcFunction)==0)
		ThrowSAXException("Boundary condition function of position is empty");

	Expression *function = Expression::CreateExpression(bcFunction,"Boundary condition function of position is not valid");
	
	unordered_map<string,double> vars;
	vars["x"] = nd[nodeNum]->x;
	vars["y"] = nd[nodeNum]->y;
	vars["r"] = nd[nodeNum]->x;
	vars["z"] = nd[nodeNum]->y;
	try
	{	bcValue = function->EvaluateFunction(vars);
	}
	catch(...)
	{	ThrowSAXException("Boundary condition function of position is not valid");
	}
	delete function;
}

// return pointer to the bcValue - a double
double *FEABoundaryCondition::GetValuePtr(void) { return &bcValue; }

