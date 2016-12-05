/********************************************************************************
    ConstraintController.cpp
    NairnFEA
    
    Created by John Nairn on 2/6/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_FEA/ConstraintController.hpp"
#include "Boundary_Conditions/Constraint.hpp"

ConstraintController *constraintCtrl=NULL;

#pragma mark Methods

// add object to linked list, but only if unique
bool ConstraintController::AddObject(LinkedObject *obj)
{
	numConstraints++;
	((Constraint *)obj)->SetLambdaNum(numConstraints);
	
	// add new one
	return ParseController::AddObject(obj);
}

