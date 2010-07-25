/********************************************************************************
	InitializationTask.cpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	Input Variables
		none
 
	Output Variables
		mpm[]->pFext
		thermal.isoRamp
********************************************************************************/

#include "NairnMPM_Class/InitializationTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Exceptions/MPMWarnings.hpp"
#include "Boundary_Conditions/MatPtLoadBC.hpp"
#include "Cracks/CrackNode.hpp"
#include "Global_Quantities/ThermalRamp.hpp"

#pragma mark CONSTRUCTORS

InitializationTask::InitializationTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get mass matrix, find dimensionless particle locations,
//	and find grid momenta
void InitializationTask::Execute(void)
{
#ifdef _PROFILE_TASKS_
	double beginTime=fmobj->CPUTime();
#endif
	
	int i;
	
	// Zero Mass Matrix and vectors
	warnings.BeginStep();
	
	// zero all nodal variabels
    for(i=1;i<=nnodes;i++)
        nd[i]->InitializeForTimeStep();
    
    // Update forces applied to particles
	MatPtLoadBC::SetParticleFext(mtime);
	
	// remove contact conditions
	CrackNode::RemoveCrackNodes();
	
    // turn off isothermal ramp when done and ramp step initialization
	thermal.CheckDone(mtime);
	
#ifdef _PROFILE_TASKS_
	totalTaskTime+=fmobj->CPUTime()-beginTime;
#endif
}	
