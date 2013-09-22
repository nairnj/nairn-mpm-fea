/********************************************************************************
	MoveCracksTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	First move crack surface in their local velocity field.
 
	Next move crack plane by one of two methods
		1. If !contact.GetMoveOnlySurfaces(), then get CM velocities for particles
			and move crack plane in the CM field
		2. Otherwise move to mean for above and below plane positions
		After moving, optionally verify no surface has crossed a crack plane
 
	Undate crack tractions
 
	Input Variables
		?

	Output Variables
		?
********************************************************************************/

#include "NairnMPM_Class/MoveCracksTask.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark CONSTRUCTORS

MoveCracksTask::MoveCracksTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Run this task to move cracks and update tactions
// throws CommonException() if crack surface or position has left the grid
void MoveCracksTask::Execute(void)
{
	// move crack surface in their local velocity fields
	CrackHeader *nextCrack=firstCrack;
	while(nextCrack!=NULL)
	{   if(!nextCrack->MoveCrack(ABOVE_CRACK))
		{	char errMsg[100];
			sprintf(errMsg,"Crack No. %d surface (above) has left the grid.",nextCrack->GetNumber());
			throw CommonException(errMsg,"NairnMPM::MPMStep");
		}	
		if(!nextCrack->MoveCrack(BELOW_CRACK))
		{	char errMsg[100];
			sprintf(errMsg,"Crack No. %d surface (below) has left the grid.",nextCrack->GetNumber());
			throw CommonException(errMsg,"NairnMPM::MPMStep");
		}
		nextCrack=(CrackHeader *)nextCrack->GetNextObject();
	}
	
	// if moving crack plane in cm velocity field, get those velocities now
	if(!contact.GetMoveOnlySurfaces()) NodalPoint::GetGridCMVelocitiesTask8();
	
	// Move crack plane by one of two methods. When moving only surface, the plane will move to average
	// of the two surfaces
	// After moving, crack plane is check for crossing, if that option is activated
	nextCrack=firstCrack;
	while(nextCrack!=NULL)
	{   if(!nextCrack->MoveCrack())
		{	char errMsg[100];
			sprintf(errMsg,"Crack No. %d position or surface has left the grid.",nextCrack->GetNumber());
			throw CommonException(errMsg,"NairnMPM::MPMStep");
		}
		nextCrack=(CrackHeader *)nextCrack->GetNextObject();
	}
	
	// update crack tractions
	if(fmobj->hasTractionCracks)
	{	nextCrack=firstCrack;
		while(nextCrack!=NULL)
		{	nextCrack->UpdateCrackTractions();
			nextCrack=(CrackHeader *)nextCrack->GetNextObject();
		}
	}
	
}
	
