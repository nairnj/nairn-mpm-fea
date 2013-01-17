/********************************************************************************
	ResetElementsTask.cpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	Check if any particles have left their element and if yes, find the new
		element. If particle reaches edge of the grid either stop the analysis
		or push it back and try to continue.
********************************************************************************/

#include "NairnMPM_Class/ResetElementsTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Exceptions/MPMTermination.hpp"
#include "Exceptions/MPMWarnings.hpp"

#pragma mark CONSTRUCTORS

ResetElementsTask::ResetElementsTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// See if any particles have changed elements
// Stop if off the grid
void ResetElementsTask::Execute(void)
{
#ifdef _PROFILE_TASKS_
	double beginTime=fmobj->CPUTime();
#endif
	int p;
	
    for(p=0;p<nmpms;p++)
    {	if(!ResetElement(mpm[p]))
		{	if(warnings.Issue(fmobj->warnParticleLeftGrid,-1)==REACHED_MAX_WARNINGS)
            {   mpm[p]->Describe();
				char errMsg[100];
			    sprintf(errMsg,"Particle No. %d left the grid\n  (plot x displacement to see it).",p+1);
				mpm[p]->origpos.x=-1.e6;
				throw MPMTermination(errMsg,"ResetElementsTask::Execute");
			}
		
			// bring back to the previous element
			ReturnToElement(mpm[p]);
		}
    }
	
#ifdef _PROFILE_TASKS_
	totalTaskTime+=fmobj->CPUTime()-beginTime;
#endif
}

// Find element for particle. Return FALSE if left
//		the grid or for GIMP moved to an edge element
//		or for axisymmetric if as R<0
// Also used during initialization to set particle's initial element
int ResetElementsTask::ResetElement(MPMBase *mpt)
{
    // check current element
    if(theElements[mpt->ElemID()]->PtInElement(mpt->pos)) return TRUE;
	
	// check neighbors if possible
	int i=0,j,elemNeighbors[27];
	theElements[mpt->ElemID()]->GetListOfNeighbors(elemNeighbors);
	while(elemNeighbors[i]!=0)
	{	j=elemNeighbors[i]-1;
    	if(theElements[j]->PtInElement(mpt->pos))
		{	if(theElements[j]->OnTheEdge()) return FALSE;
			if(fmobj->IsAxisymmetric() && mpt->pos.x<0.) return FALSE;
			mpt->ChangeElemID(j);
			return TRUE;
		}
		i++;
    }
    
    // if still not found, check all elements
    for(i=0;i<nelems;i++)
    {	if(theElements[i]->PtInElement(mpt->pos))
		{	if(theElements[i]->OnTheEdge()) return FALSE;
			if(fmobj->IsAxisymmetric() && mpt->pos.x<0.) return FALSE;
			mpt->ChangeElemID(i);
			return TRUE;
		}
    }
    return FALSE;
}
	
// Find element for particle. Return FALSE if left
// the grid or for GIMP moved to an edge element
void ResetElementsTask::ReturnToElement(MPMBase *mpt)
{
	Vector outside=mpt->pos;
    int elemID=mpt->ElemID();
	Vector inside,middle,origin;
	int pass;
	
	// try to retrace position, if fails, take element centroid
	inside.x=outside.x-timestep*mpt->vel.x;
	inside.y=outside.y-timestep*mpt->vel.y;
	inside.z=outside.z-timestep*mpt->vel.z;
	if(!theElements[elemID]->PtInElement(inside))
		theElements[elemID]->GetXYZCentroid(&inside);
	origin=inside;
	
	// bisect 10 times
	for(pass=1;pass<=10;pass++)
	{	middle.x=(outside.x+inside.x)/2.;
		middle.y=(outside.y+inside.y)/2.;
		middle.z=(outside.z+inside.z)/2.;
		
		if(theElements[elemID]->PtInElement(middle))
			inside=middle;
		else
			outside=middle;
	}
	
	// move to inside
	mpt->SetPosition(&inside);
	
	// change velocity for movement form starting position to new edge position (but seems to not be good idea)
	//mpt->vel.x=(inside.x-origin.x)/timestep;
	//mpt->vel.y=(inside.y-origin.y)/timestep;
	//mpt->vel.z=(inside.z-origin.z)/timestep;
}



