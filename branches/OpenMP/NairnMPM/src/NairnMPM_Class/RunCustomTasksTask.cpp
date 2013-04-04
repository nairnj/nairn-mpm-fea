/********************************************************************************
	RunCustomTasksTask.cpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	Run all custom tasks, but first run loop that allows then to extrapolate
		to the grid
********************************************************************************/

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
	// zero thisfunction in case in 2D analysis
	for(int i=0;i<MaxShapeNds;i++) zDeriv[i]=0.;
}

#pragma mark REQUIRED METHODS

// Run all custom tasks
void RunCustomTasksTask::Execute(void)
{
	int i,p,iel,matfld,numnds,nds[MaxShapeNds];
	MaterialBase *matID;
	double wt,fn[MaxShapeNds],xDeriv[MaxShapeNds],yDeriv[MaxShapeNds];
	short vfld,isRigid;
	
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
	*/
    if(needExtrapolations)
    {	// call each task for initialization prior to extrapolations
    	nextTask=theTasks;
        while(nextTask!=NULL)
            nextTask=nextTask->BeginExtrapolations();
        
        // particle loop
        for(p=0;p<nmpms;p++)
        {   // Load element coordinates
			matID=theMaterials[mpm[p]->MatID()];
			if(matID->RigidBC()) continue;			// skip boundary condition particle always
			isRigid=matID->Rigid();					// if TRUE, will be rigid contact particle
			matfld=matID->GetField();
			
            // find shape functions and derviatives
			iel=mpm[p]->ElemID();
			theElements[iel]->GetShapeGradients(&numnds,fn,nds,mpm[p]->GetNcpos(),xDeriv,yDeriv,zDeriv,mpm[p]);
            
			// Add particle property to each node in the element
            for(i=1;i<=numnds;i++)
            {   // global mass matrix
				vfld=(short)mpm[p]->vfld[i];				// velocity field to use
                wt=fn[i]*mpm[p]->mp;
                
                // possible extrapolation to the nodes
                nextTask=theTasks;
                while(nextTask!=NULL)
                    nextTask=nextTask->NodalExtrapolation(nd[nds[i]],mpm[p],vfld,matfld,wt,isRigid);
				
                // possible extrapolation to the particle (but currently not used by any custom task)
                /*
                nextTask=theTasks;
                while(nextTask!=NULL)
                    nextTask=nextTask->ParticleCalculation(nd[nds[i]],mpm[p],vfld,matfld,fn[i],xDeriv[i],yDeriv[i],zDeriv[i],isRigid);
                */
            }
            
            // possible single calculations for each particle
            nextTask=theTasks;
            while(nextTask!=NULL)
                nextTask=nextTask->ParticleExtrapolation(mpm[p],isRigid);
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
    nextTask=theTasks;
    while(nextTask!=NULL)
    	nextTask=nextTask->FinishForStep();
	
}
