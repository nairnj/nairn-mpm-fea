/********************************************************************************
	ResetElementsTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	Check if any particles have left their element and if yes, find the new
		element. If particle reaches edge of the grid either stop the analysis
		or push it back and try to continue.
 
	Also update alpha for feedback damping. Wanted to have this at the end
		of time step fro both USF and USAVG. This update used to be in the
		UpdateParticlesTask, but that had problems when using grid kinetic
		energy of cohesive zones. Here is betters, especially when using
		USAVG mode.
********************************************************************************/

#include "NairnMPM_Class/ResetElementsTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Exceptions/MPMWarnings.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "Patches/GridPatch.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark CONSTRUCTORS

ResetElementsTask::ResetElementsTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// See if any particles have changed elements
// Stop if off the grid
void ResetElementsTask::Execute(void)
{
	CommonException *resetErr = NULL;
	
	// update feedback damping now if needed
	bodyFrc.UpdateAlpha(timestep,mtime);

	// This block should be made parallel
	// But when do so, need method to move particle  between patches while keeping threads independent
	int totalPatches = fmobj->GetTotalNumberOfPatches();
	
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
					{	if(warnings.Issue(fmobj->warnParticleLeftGrid,-1)==REACHED_MAX_WARNINGS)
						{	// print message and log error
							if(resetErr==NULL)
							{	mptr->Describe();
								char errMsg[100];
								sprintf(errMsg,"Too many particles have left the grid\n  (plot x displacement to see last one).");
								mptr->origpos.x=-1.e6;
								resetErr = new CommonException(errMsg,"ResetElementsTask::Execute");
								throw *resetErr;
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
				
				// next material point and update previous particle
				prevMptr = mptr;
				mptr = (MPMBase *)mptr->GetNextObject();
			}
		}
	}
	
	// if error occurred then throw it
	if(resetErr) throw *resetErr;
}

// Find element for particle. Return FALSE if left
//		the grid or for GIMP moved to an edge element
//		or for axisymmetric if as R<0
// Also used during initialization to set particle's initial element
int ResetElementsTask::ResetElement(MPMBase *mpt)
{
    // check current element
    if(theElements[mpt->ElemID()]->PtInElement(mpt->pos))
	{	// it has not changed elements
		return SAME_ELEMENT;
	}
    
/*
	// check neighbors if possible
	int i=0,j,elemNeighbors[27];
	theElements[mpt->ElemID()]->GetListOfNeighbors(elemNeighbors);
	while(elemNeighbors[i]!=0)
	{	j=elemNeighbors[i]-1;
    	if(theElements[j]->PtInElement(mpt->pos))
		{	if(theElements[j]->OnTheEdge()) return LEFT_GRID;
			if(fmobj->IsAxisymmetric() && mpt->pos.x<0.) return LEFT_GRID;
			mpt->ChangeElemID(j);
			return NEW_ELEMENT;
		}
		i++;
    }
    
    // if still not found, check all elements
    for(i=0;i<nelems;i++)
    {	if(theElements[i]->PtInElement(mpt->pos))
		{	if(theElements[i]->OnTheEdge()) return LEFT_GRID;
			if(fmobj->IsAxisymmetric() && mpt->pos.x<0.) return LEFT_GRID;
			mpt->ChangeElemID(i);
			return NEW_ELEMENT;
		}
    }
 
    return LEFT_GRID;
*/
	
    if(mpmgrid.IsStructuredGrid())
    {   try
        {   // calculate from coordinates
            int j = mpmgrid.FindElementFromPoint(&mpt->pos)-1;
            if(theElements[j]->OnTheEdge()) return LEFT_GRID;
			if(fmobj->IsAxisymmetric() && mpt->pos.x<0.) return LEFT_GRID;
            mpt->ChangeElemID(j);
            return NEW_ELEMENT;
        }
        catch(...)
        {   return LEFT_GRID;
        }
    }
	
    // for unstructured grid, have to check all elements
    for(int i=0;i<nelems;i++)
    {	if(theElements[i]->PtInElement(mpt->pos))
        {	if(theElements[i]->OnTheEdge()) return LEFT_GRID;
            if(fmobj->IsAxisymmetric() && mpt->pos.x<0.) return LEFT_GRID;
            mpt->ChangeElemID(i);
            return NEW_ELEMENT;
        }
    }
    
    return LEFT_GRID;
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



