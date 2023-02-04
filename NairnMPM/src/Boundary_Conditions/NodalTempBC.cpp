/********************************************************************************
    NodalTempBC.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Oct 18 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Boundary_Conditions/NodalTempBC.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"

// Nodal temperature BC global
NodalTempBC *firstTempBC=NULL;
NodalTempBC *lastTempBC=NULL;
NodalTempBC *firstRigidTempBC=NULL;
NodalTempBC *reuseRigidTempBC=NULL;

#pragma mark NodalTempBC: Constructors and Destructors

NodalTempBC::NodalTempBC(int num,int setStyle,double temperature,double argTime)
		: NodalValueBC(num,setStyle,temperature,argTime)
{
	temperatureNoBC = NULL;
}

#pragma mark NodalTempBC: Methods

// save nodal temperature and zero it
// throws std::bad_alloc
NodalTempBC *NodalTempBC::CopyNodalValue(NodalPoint *nd,TransportField *gTrans)
{
	// create vector to hold options
	if(temperatureNoBC==NULL)
	{
		temperatureNoBC = new double[1];
	}
	
	// copy global one first
	temperatureNoBC[0] = nd->gCond.gTValue;

	// return next task
	return (NodalTempBC *)GetNextObject();
}

// restore nodal temperature and get initial force to cancel no-BC result
NodalTempBC *NodalTempBC::PasteNodalValue(NodalPoint *nd,TransportField *gTrans)
{
	// paste global temperature
	nd->gCond.gTValue = temperatureNoBC[0];

	// return next task
	return (NodalTempBC *)GetNextObject();
}

// when getting total heat reaction, add qreaction to input double
// if matchID==0 include it, otherwise include only if matchID equals bcID
NodalTempBC *NodalTempBC::AddHeatReaction(double *totalReaction,int matchID)
{	if(bcID==matchID || matchID==0)
		*totalReaction += qreaction;
	return (NodalTempBC *)GetNextObject();
}

#pragma mark ACCESSORS

// get set direction
int NodalTempBC::GetSetDirection(void) const { return TEMP_DIRECTION; }

/**********************************************************
	Sum all reaction heat forces for all temperature BCs
	If ID is not zero, only includes those with matching ID
	Result is in nJ/sec when in Legacy units
*/
double NodalTempBC::TotalHeatReaction(int matchID)
{
	double reactionTotal = 0.;
    NodalTempBC *nextBC=firstTempBC;
    
    // Add each value
    while(nextBC!=NULL)
		nextBC=nextBC->AddHeatReaction(&reactionTotal,matchID);
	
	return reactionTotal;
}

