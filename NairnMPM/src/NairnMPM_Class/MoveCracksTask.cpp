/********************************************************************************
	MoveCracksTask.cpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	Move cracks and crack surfaces.
********************************************************************************/

#include "NairnMPM_Class/MoveCracksTask.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Exceptions/MPMTermination.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"

#pragma mark CONSTRUCTORS

MoveCracksTask::MoveCracksTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Run all custom tasks
void MoveCracksTask::Execute(void)
{
#ifdef _PROFILE_TASKS_
	double beginTime=fmobj->CPUTime();
#endif

	// prepare multimaterial fields for moving cracks
	CrackHeader *nextCrack=firstCrack;
	while(nextCrack!=NULL)
	{   if(!nextCrack->MoveCrack(ABOVE_CRACK))
		{	char errMsg[100];
			sprintf(errMsg,"Crack No. %d surface (above) has left the grid.",nextCrack->GetNumber());
			throw MPMTermination(errMsg,"NairnMPM::MPMStep");
		}	
		if(!nextCrack->MoveCrack(BELOW_CRACK))
		{	char errMsg[100];
			sprintf(errMsg,"Crack No. %d surface (below) has left the grid.",nextCrack->GetNumber());
			throw MPMTermination(errMsg,"NairnMPM::MPMStep");
		}
		nextCrack=(CrackHeader *)nextCrack->GetNextObject();
	}
	
	if(!contact.GetMoveOnlySurfaces()) NodalPoint::GetGridCMVelocitiesTask8();
	nextCrack=firstCrack;
	while(nextCrack!=NULL)
	{   if(!nextCrack->MoveCrack())
		{	char errMsg[100];
			sprintf(errMsg,"Crack No. %d position or surface has left the grid.",nextCrack->GetNumber());
			throw MPMTermination(errMsg,"NairnMPM::MPMStep");
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
	
#ifdef _PROFILE_TASKS_
	totalTaskTime+=fmobj->CPUTime()-beginTime;
#endif
}
	
