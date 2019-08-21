/********************************************************************************
	ResetElementsTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	The tasks are:
	--------------
	* Update alpha for feedback damping. Wanted to have this at the end
      of time step for both USF and USAVG. This update used to be in the
      UpdateParticlesTask, but that had problems when using grid kinetic
      energy of cohesive zones. Here is betters, especially when using
      USAVG mode.
	* Check if any particles have left their element and if yes, find the new
      element.
	* If particle reaches edge of the grid either stop the analysis
	  or push it back and try to continue.
 
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/ResetElementsTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Exceptions/MPMWarnings.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "Patches/GridPatch.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Exceptions/CommonException.hpp"

// uncomment to make this task parallel
#define PARALLEL_RESET

#pragma mark CONSTRUCTORS

ResetElementsTask::ResetElementsTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// See if any particles have changed elements
// Stop if off the grid
// throws CommonException()
void ResetElementsTask::Execute(int taskOption)
{
	// update feedback damping now if needed
	bodyFrc.UpdateAlpha(timestep,mtime);
	
	// how many patches?
	int totalPatches = fmobj->GetTotalNumberOfPatches();
	
#ifdef PARALLEL_RESET
	// initialize error
	CommonException *resetErr = NULL;

	// parallel over patches
#pragma omp parallel
	{
        // thread for patch pn
		int pn = GetPatchNumber();
        
		try
		{	// resetting all element types
			for(int block=FIRST_NONRIGID;block<=FIRST_RIGID_BC;block++)
			{	// get first material point in this block
				MPMBase *mptr = patches[pn]->GetFirstBlockPointer(block);
				MPMBase *prevMptr = NULL;		// previous one of this type in current patch
				while(mptr!=NULL)
				{	int status = ResetElement(mptr);
					
					if(status==LEFT_GRID)
					{	// particle has left the grid
						mptr->IncrementElementCrossings();
						
						// enter warning only if this particle did not leave the grid before
						if(!mptr->HasLeftTheGridBefore())
						{	int result = warnings.Issue(fmobj->warnParticleLeftGrid,-1);
							if(result==REACHED_MAX_WARNINGS || result==GAVE_WARNING)
							{
#pragma omp critical (output)
								{	mptr->Describe();
								}
								// abort if needed
								if(result==REACHED_MAX_WARNINGS)
								{	char errMsg[100];
									sprintf(errMsg,"Too many particles have left the grid\n  (plot x displacement to see last one).");
									mptr->origpos.x=-1.e6;
									throw CommonException(errMsg,"ResetElementsTask::Execute");
								}
							}
							
							// set this particle has left the grid once
							mptr->SetHasLeftTheGridBefore(TRUE);
						}
						
						// bring back to the previous element
						ReturnToElement(mptr);
					}
					
					else if(status==NEW_ELEMENT && totalPatches>1)
					{	// did it also move to a new patch?
						int newpn = mpmgrid.GetPatchForElement(mptr->ElemID());
						if(pn != newpn)
						{	if(!patches[pn]->AddMovingParticle(mptr,patches[newpn],prevMptr))
							{	throw CommonException("Out of memory storing data for particle changing patches","ResetElementsTask::Execute");
							}
						}
					}
					
					else if(status==LEFT_GRID_NAN)
					{
#pragma omp critical (output)
						{	cout << "# Particle has left the grid and position is nan" << endl;
							mptr->Describe();
						}
						throw CommonException("Particle has left the grid and position is nan","ResetElementsTask::Execute");
					}
					
					// next material point and update previous particle
					prevMptr = mptr;
					mptr = (MPMBase *)mptr->GetNextObject();
				}
			}
		}
		catch(CommonException& err)
        {   if(resetErr==NULL)
			{
#pragma omp critical (error)
				resetErr = new CommonException(err);
			}
		}
		catch(CommonException* err)
		{	if(resetErr==NULL)
			{
#pragma omp critical (error)
				resetErr = new CommonException(*err);
			}
		}
		catch(std::bad_alloc&)
        {   if(resetErr==NULL)
			{
#pragma omp critical (error)
				resetErr = new CommonException("Memory error","ResetElementsTask::Execute");
			}
		}
		catch(...)
        {   if(resetErr==NULL)
			{
#pragma omp critical (error)
				resetErr = new CommonException("Unexepected error","ResetElementsTask::Execute");
			}
		}
	}

	// throw now if was an error
	if(resetErr!=NULL) throw *resetErr;
    
	// reduction phase moves the particles
	for(int pn=0;pn<totalPatches;pn++)
		patches[pn]->MoveParticlesToNewPatches();
	
#else
	
	int status;
	MPMBase *mptr,*prevMptr,*nextMptr;
	for(int pn=0;pn<totalPatches;pn++)
	{	for(int block=FIRST_NONRIGID;block<=FIRST_RIGID_BC;block++)
		{	// get first material point in this block
			mptr = patches[pn]->GetFirstBlockPointer(block);
			prevMptr = NULL;		// previous one of this type in current patch
			while(mptr!=NULL)
			{	status = ResetElement(mptr);
				
				if(status==LEFT_GRID)
				{	// particle has left the grid
					mptr->IncrementElementCrossings();
				
					// enter warning only if this particle did not leave the grid before
					if(!mptr->HasLeftTheGridBefore())
					{	int result = warnings.Issue(fmobj->warnParticleLeftGrid,-1);
						if(result==REACHED_MAX_WARNINGS || result==GAVE_WARNING)
						{	mptr->Describe();
							// abort if needed
							if(result==REACHED_MAX_WARNINGS)
							{	char errMsg[100];
								sprintf(errMsg,"Too many particles have left the grid\n  (plot x displacement to see last one).");
								mptr->origpos.x=-1.e6;
								throw CommonException(errMsg,"ResetElementsTask::Execute");
							}
						}
						
						// set this particle has left the grid once
						mptr->SetHasLeftTheGridBefore(TRUE);
					}
				
					// bring back to the previous element
					ReturnToElement(mptr);
				}
				
				else if(status==NEW_ELEMENT && totalPatches>1)
				{	int newpn = mpmgrid.GetPatchForElement(mptr->ElemID());
					if(pn != newpn)
					{	// next material point read before move this particle
						nextMptr = (MPMBase *)mptr->GetNextObject();
						
						// move particle mptr
						patches[pn]->RemoveParticleAfter(mptr,prevMptr);
						patches[newpn]->AddParticle(mptr);
						
						// next material point is now after the prevMptr, which stays the same, which may be NULL
						mptr = nextMptr;
						continue;
					}
				}
				
				else if(status==LEFT_GRID_NAN)
				{	cout << "# Particle has left the grid and position is nan" << endl;
					mptr->Describe();
					throw CommonException("Particle has left the grid and position is nan","ResetElementsTask::Execute");
				}
				
				// next material point and update previous particle
				prevMptr = mptr;
				mptr = (MPMBase *)mptr->GetNextObject();
			}
		}
	}
#endif
}

// Find element for particle. Return FALSE if left
//		the grid or for GIMP moved to an edge element
//		or for axisymmetric if as R<0
// Also used during initialization to set particle's initial element
int ResetElementsTask::ResetElement(MPMBase *mpt)
{
	// check NAN
	if(mpt->pos.x!=mpt->pos.x || mpt->pos.y!=mpt->pos.y || mpt->pos.z!=mpt->pos.z)
	{	return LEFT_GRID_NAN;
	}
	
	// check current element
	if(theElements[mpt->ElemID()]->PtInElement(mpt->pos))
	{	// it has not changed elements
		return SAME_ELEMENT;
	}
    
	// calculate from coordinates
	try
	{   int j = mpmgrid.FindElementFromPoint(&mpt->pos,mpt)-1;		// elem ID (0 based)
		if(theElements[j]->OnTheEdge()) return LEFT_GRID;
		if(fmobj->IsAxisymmetric() && mpt->pos.x<=0.) return LEFT_GRID;
		mpt->ChangeElemID(j,!mpmgrid.IsStructuredEqualElementsGrid());
		return NEW_ELEMENT;
	}
	catch(...)
	{   return LEFT_GRID;
	}
}
	
// Push particle back to its previous element
// the grid or for GIMP moved to an edge element
void ResetElementsTask::ReturnToElement(MPMBase *mpt)
{
	Vector outside = mpt->pos;
    int elemID = mpt->ElemID();
	Vector inside,middle;
	int pass;
	
	// try to retrace position, if fails, take element centroid
	inside.x = outside.x - timestep*mpt->vel.x;
	inside.y = outside.y - timestep*mpt->vel.y;
	inside.z = outside.z - timestep*mpt->vel.z;
	if(!theElements[elemID]->PtInElement(inside))
		theElements[elemID]->GetXYZCentroid(&inside);
	
	// bisect 10 times until inside, if fails, uses above inside
	for(pass=1;pass<=10;pass++)
	{	middle.x = (outside.x+inside.x)/2.;
		middle.y = (outside.y+inside.y)/2.;
		middle.z = (outside.z+inside.z)/2.;
		
		if(theElements[elemID]->PtInElement(middle))
			inside=middle;
		else
			outside=middle;
	}
	
	// move to inside
	mpt->SetPosition(&inside);
	
	// change velocity for movement from starting position to new edge position (but seems to not be good idea)
	//mpt->vel.x=(inside.x-origin.x)/timestep;
	//mpt->vel.y=(inside.y-origin.y)/timestep;
	//mpt->vel.z=(inside.z-origin.z)/timestep;
}



