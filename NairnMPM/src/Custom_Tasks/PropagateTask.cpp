/********************************************************************************
    PropagateTask.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Mon Aug 18 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Custom_Tasks/PropagateTask.hpp"
#include "Custom_Tasks/CalcJKTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Cracks/CrackSegment.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "System/ArchiveData.hpp"
#include "System/UnitsController.hpp"

// globals
PropagateTask *propagateTask=NULL;
double PropagateTask::cellsPerPropagationStep=0.5;

#pragma mark Constructors and Destructors

// Constructors
PropagateTask::PropagateTask()
{
    propagateTask = this;
    arrested = false;
}

// Return name of this task
const char *PropagateTask::TaskName(void) { return "Crack Propagation Task"; }

#pragma mark GENERIC TASK METHODS

// at beginning of analysis
CustomTask *PropagateTask::Initialize(void)
{
    nextPropTime=propTime;
    cout << "Crack propagation activated." << endl;
	cout << "   Crack growth by steps = " << cellsPerPropagationStep <<"*(cell size)" << endl;
    return nextTask;
}

// called when MPM step is getting ready to do custom tasks
// This tasks only exists if propagation is turned on
CustomTask *PropagateTask::PrepareForStep(bool &needExtraps)
{
    theResult=NOGROWTH;
    
    if(mtime+timestep>=nextPropTime && !arrested)
    {	// ready to do crack propagation
    	doPropCalcs = true;
        nextPropTime += propTime;
        
        // Make sure J and K are available if needed for current criterion
		CrackHeader *nextCrack=firstCrack;
		while(nextCrack!=NULL)
		{	if(theJKTask!=NULL) theJKTask->ScheduleJK(nextCrack->CriterionNeeds());
			nextCrack=(CrackHeader *)nextCrack->GetNextObject();
		}
    }
    else
        doPropCalcs = false;
        
    return nextTask;
}

// Calculate J and K at crack tips
CustomTask *PropagateTask::StepCalculation(void)
{
    // if not needed, just exit
    if(!doPropCalcs) return nextTask;
	
    CrackHeader *nextCrack;
    CrackSegment *crkTip;
    double cSize;
    int i,inMat;
    Vector tipDir,growTo,grow;
    int tipElem;
	char isAlt[10];
    
    // loop over cracks
    nextCrack=firstCrack;
    while(nextCrack!=NULL)
    {	// each crack tip
        for(i=START_OF_CRACK;i<=END_OF_CRACK;i++)
        {   // find crack tip and direction
			nextCrack->CrackTipAndDirection(i,&crkTip,tipDir);
            
            // crack propagation
            inMat=crkTip->tipMatnum;
            if(inMat>0)
            {	// see if it grows
				int shouldGo=theMaterials[inMat-1]->ShouldPropagate(crkTip,tipDir,nextCrack,fmobj->np,0);
				isAlt[0] = 0;
				if(shouldGo==GROWNOW)
				{	nextCrack->SetAllowAlternate(i,false);
				}
				else if(nextCrack->GetAllowAlternate(i))
				{	shouldGo=theMaterials[inMat-1]->ShouldPropagate(crkTip,tipDir,nextCrack,fmobj->np,1);
					if(shouldGo==GROWNOW) strcpy(isAlt," (alt)");
				}
				if(shouldGo==GROWNOW)
                {   theResult=GROWNOW;
                    tipElem=crkTip->planeElemID();
                    if(fabs(tipDir.x)>fabs(tipDir.y))
                        cSize=theElements[tipElem]->xmax-theElements[tipElem]->xmin;
                    else
                        cSize=theElements[tipElem]->ymax-theElements[tipElem]->ymin;
                    grow.x=cellsPerPropagationStep*cSize*tipDir.x;
                    grow.y=cellsPerPropagationStep*cSize*tipDir.y;
					
					// adjust if crossing another crack - adjusts grow is need and returns resulting relative change (in p)
					double p = nextCrack->AdjustGrowForCrossing(&grow,crkTip,cSize,&tipDir);
					
					// crack number and tip
					archiver->IncrementPropagationCounter();
                    cout << "# propagation" << isAlt << " crack-tip " << nextCrack->GetNumber() << "-" << i;
					
					// summarize
					cout << " at t=" << mtime*UnitsController::Scaling(1.e3)
						<< " with J=Jtip+Jzone : " << crkTip->Jint.z*UnitsController::Scaling(1.e-3)
						<< " = " << crkTip->Jint.x*UnitsController::Scaling(1.e-3)
						<< " + " << (crkTip->Jint.z-crkTip->Jint.x)*UnitsController::Scaling(1.e-3) << endl;
                    
                    // if jump is .7 or more cells, divide propagation into multiple segments
                    int iseg,numSegs = 1;
                    if(p*cellsPerPropagationStep>.7) numSegs= (int)(2*(p*cellsPerPropagationStep+.25));
                    CrackSegment *newCrkTip = NULL;
                    for(iseg=1;iseg<=numSegs;iseg++)
                    {   growTo.x=crkTip->cp.x+(double)iseg*grow.x/(double)numSegs;
                        growTo.y=crkTip->cp.y+(double)iseg*grow.y/(double)numSegs;
						
						// if crack has traction propagation, use it, otherwise use material point
						int tractionID = nextCrack->GetTractionPropID();
						if(tractionID<=0)
						{	tractionID = isAlt[0]==0 ? theMaterials[inMat-1]->tractionMat[0] :
														theMaterials[inMat-1]->tractionMat[1] ;
						}
						newCrkTip=nextCrack->Propagate(growTo,(int)i,tractionID);
					}
                    crkTip = newCrkTip;
                    
					if(crkTip!=NULL)
					{	// check if crack speed is being controlled
						if(theMaterials[inMat-1]->ControlCrackSpeed(crkTip,propTime))
							nextPropTime=mtime+propTime;
						
						// crack tip heating (if activated)
						if(ConductionTask::active)
							conduction->StartCrackTipHeating(crkTip,grow,nextCrack->GetThickness());
					}
                }
            }
        }

        // next crack
        nextCrack=(CrackHeader *)nextCrack->GetNextObject();
    }
    
    return nextTask;
}

#pragma mark PropagateTask METHODS

// turn crack propagation task on or off
void PropagateTask::ArrestGrowth(bool newArrest) { arrested=newArrest; }
bool PropagateTask::Arrested(void) { return arrested; }

