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
	CommonException *mcErr = NULL;
    
    int numProcsPerCrack = (double)numberOfCracks/(double)(fmobj->GetNumberOfProcessors());
    
    // uncomment this line to prevent parallelization of this task
    //numProcsPerCrack = 0;
    
	// move crack surface in their local velocity fields
#pragma omp parallel for if(numProcsPerCrack>1)
    for(int cn=0;cn<numberOfCracks;cn++)
    {   try
        {   if(!crackList[cn]->MoveCrack(ABOVE_CRACK))
            {	char errMsg[100];
                sprintf(errMsg,"Crack No. %d surface (above) has left the grid.",cn+1);
                throw CommonException(errMsg,"MoveCracksTask::Execute");
            }
            if(!crackList[cn]->MoveCrack(BELOW_CRACK))
            {	char errMsg[100];
                sprintf(errMsg,"Crack No. %d surface (below) has left the grid.",cn+1);
                throw CommonException(errMsg,"MoveCracksTask::Execute");
            }
        }
		catch(CommonException err)
		{	if(mcErr==NULL)
            {
#pragma omp critical
                mcErr = new CommonException(err);
            }
		}
	}
    
	// throw any errors
	if(mcErr!=NULL) throw *mcErr;
    
	// if moving crack plane in ctr mass velocity field, get those velocities now
	if(!contact.GetMoveOnlySurfaces()) NodalPoint::GetGridCMVelocitiesTask8();
	
	// Move crack plane by one of two methods. When moving only surface, the plane will move to average
	// of the two surfaces found above, otherwise moves in the just-calculated center or mass field
	// After moving, crack plane is checked for crossing, if that option is activated
#pragma omp parallel for if(numProcsPerCrack>1)
    for(int cn=0;cn<numberOfCracks;cn++)
    {   try
        {   if(!crackList[cn]->MoveCrack())
            {	char errMsg[100];
                sprintf(errMsg,"Crack No. %d position or surface has left the grid.",cn+1);
                throw CommonException(errMsg,"MoveCracksTask::Execute");
            }
        }
		catch(CommonException err)
		{	if(mcErr==NULL)
            {
#pragma omp critical
                mcErr = new CommonException(err);
            }
		}
	}
	
	// throw any errors
	if(mcErr!=NULL) throw *mcErr;
    
	// update crack tractions
	if(fmobj->hasTractionCracks)
	{
#pragma omp parallel for if(numProcsPerCrack>1)
        for(int cn=0;cn<numberOfCracks;cn++)
		{	crackList[cn]->UpdateCrackTractions();
		}
	}
    
}
	
