/********************************************************************************
    NodalConcBC.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Apr 1 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Boundary_Conditions/NodalConcBC.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "System/UnitsController.hpp"

// Nodal concentration BC globals
NodalConcBC *firstConcBC=NULL;
NodalConcBC *lastConcBC=NULL;
NodalConcBC *firstRigidConcBC=NULL;
NodalConcBC *reuseRigidConcBC=NULL;

#pragma mark NodalConcBC: Constructors and Destructors

NodalConcBC::NodalConcBC(int num,int setStyle,double concentration,double argTime)
                    : NodalValueBC(num,setStyle,concentration,argTime)
{
}

// print it (if can be used)
BoundaryCondition *NodalConcBC::PrintBC(ostream &os)
{
	// moisture uses parent calss
	if(DiffusionTask::HasDiffusion())
		return NodalValueBC::PrintBC(os);
	
	return (BoundaryCondition *)GetNextObject();
}

#pragma mark NodalConcBC: ACCESSORS

// set value (and scale legacy MPa to Pa)
void NodalConcBC::SetBCValue(double bcvalue)
{	double rescale = DiffusionTask::RescalePotential();
	BoundaryCondition::SetBCValue(rescale*bcvalue);
}

// get set direction
int NodalConcBC::GetSetDirection(void) const { return CONC_DIRECTION; }

// return point on node to transport field
TransportField *NodalConcBC::GetTransportFieldPtr(NodalPoint *ndpt) const { return &(ndpt->gDiff); }
