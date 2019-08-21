/********************************************************************************
	RunCustomTasksTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	Run all custom tasks. The tasks area:
	------------------------------------
	* If any tasks needs extrapolations, all them to extrapolate to the grid
		- Each custom tasks uses BeginExtrapolations(), NodalExtrapolation(),
		  and Extrapolations().
	* Call StepCalculation() for each custom task
	* Call FinishForStep() for each custom task
	  (a task can delete itself if all done)
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/RunCustomTasksTask.hpp"
#include "Custom_Tasks/CustomTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"

#pragma mark CONSTRUCTORS

RunCustomTasksTask::RunCustomTasksTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Run all custom tasks
void RunCustomTasksTask::Execute(int taskOption)
{
#ifdef CONST_ARRAYS
	int ndsArray[MAX_SHAPE_NODES];
	double fn[MAX_SHAPE_NODES];
	//double xDeriv[MAX_SHAPE_NODES],yDeriv[MAX_SHAPE_NODES],zDeriv[MAX_SHAPE_NODES];
#else
	int ndsArray[maxShapeNodes];
	double fn[maxShapeNodes];
	//double xDeriv[maxShapeNodes],yDeriv[maxShapeNodes],zDeriv[maxShapeNodes];
#endif
	int matfld,numnds;
	const MaterialBase *matID;
	double fnmp;
	short vfld,isRigid;
    
	// in case 2D planar
    //for(int i=0;i<maxShapeNodes;i++) zDeriv[i] = 0.;
	
    /* Step 1: Call all tasks. The tasks can do initializations needed
			for this step. If any task needs nodal extrapolations
			they should set needNodalExtraps to TRUE. If it does
			not need them, leave it alone
	*/
    bool needExtrapolations=FALSE;
    CustomTask *nextTask=theTasks;
    while(nextTask!=NULL)
	{	bool taskNeedsExtrapolations=FALSE;
    	nextTask=nextTask->PrepareForStep(taskNeedsExtrapolations);
		// if it was set to TRUE, trasfer to to global setting
		if(taskNeedsExtrapolations) needExtrapolations=TRUE;
	}
    
    /* Step 2: Extrapolate particle info to grid if needed for
			any custom task
	   This loop is not parallel. It should be avoided except by custom tasks
			that only run periodically (such as for archiving)
	*/
    if(needExtrapolations)
    {	// call each task for initialization prior to extrapolations
    	nextTask=theTasks;
        while(nextTask!=NULL)
            nextTask=nextTask->BeginExtrapolations();
        
        // particle loop or nonrigid, rigid block, and rigid contact particles
        for(int p=0;p<nmpmsRC;p++)
		{	MPMBase *mpmptr = mpm[p];
			
           // Load element coordinates
			matID=theMaterials[mpmptr->MatID()];
			isRigid=matID->IsRigid();					// if TRUE, will be rigid contact particle
			matfld=matID->GetField();
			
            // find shape functions and derviatives
			const ElementBase *elref = theElements[mpmptr->ElemID()];
			int *nds = ndsArray;
			elref->GetShapeFunctions(fn,&nds,mpmptr);
			numnds = nds[0];
            
			// Add particle property to each node in the element
            for(int i=1;i<=numnds;i++)
            {   // global mass matrix
				vfld=(short)mpmptr->vfld[i];				// velocity field to use
                fnmp=fn[i]*mpmptr->mp;
                
                // possible extrapolation to the nodes
                nextTask=theTasks;
                while(nextTask!=NULL)
                    nextTask=nextTask->NodalExtrapolation(nd[nds[i]],mpmptr,vfld,matfld,fnmp,isRigid);
				
             }
        }
        
        // finished with extrapolations
        nextTask=theTasks;
        while(nextTask!=NULL)
            nextTask=nextTask->EndExtrapolations();
    }
	
    // Step 3: Do the custom task calculations
    nextTask=theTasks;
    while(nextTask!=NULL)
        nextTask=nextTask->StepCalculation();
	
    // Step 4: Call tasks in case any need to clean up
    CustomTask *thisTask=theTasks;
	CustomTask *prevTask=NULL;
    while(thisTask!=NULL)
	{	bool removeTask = false;
    	nextTask = thisTask->FinishForStep(removeTask);
		if(removeTask)
		{	// redirect previous task to nextTask, thus skipping thisTask
			if(prevTask==NULL)
			{	// removing the first task, so next task will now be the first
				theTasks = nextTask;
			}
			else
			{	// point previous task to next task, thus skipping this task
				((CustomTask *)prevTask)->nextTask = nextTask;
			}
			// now delete it (prevTask stays the same)
			delete thisTask;
		}
		else
		{	// this is now the previous one
			prevTask = thisTask;
		}
		// on to next task
		thisTask = nextTask;
	}
	
}
