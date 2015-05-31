/********************************************************************************
    PropagateTask.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Mon Aug 18 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

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

// globals
PropagateTask *propagateTask=NULL;
double PropagateTask::cellsPerPropagationStep=0.5;

#pragma mark Constructors and Destructors

// Constructors
PropagateTask::PropagateTask()
{
    propagateTask=this;
    arrested=FALSE;
    totalPlastic = 0.;
    totalPotential = 0.;
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
    	doPropCalcs = TRUE;
        nextPropTime += propTime;
        doEnergyBalanceCalcs = FALSE;
        
        // Make sure J and K are available if needed for current criterion
        // also see if need totalPlastic and totalPotential
		CrackHeader *nextCrack=firstCrack;
		while(nextCrack!=NULL)
		{	if(theJKTask!=NULL) theJKTask->ScheduleJK(nextCrack->CriterionNeeds(doEnergyBalanceCalcs));
			nextCrack=(CrackHeader *)nextCrack->GetNextObject();
		}
    }
    else
        doPropCalcs=FALSE;
        
    return nextTask;
}

// Calculate J and K at crack tips
CustomTask *PropagateTask::StepCalculation(void)
{
    // if not needed, just exit
    if(!doPropCalcs) return nextTask;
	
	// particle extrapolation if needed
    if(doEnergyBalanceCalcs)
    {   totalPlastic = 0.;
        totalPotential = 0.;
        for(int p=0;p<nmpmsNR;p++)
        {   MPMBase *mpnt = mpm[p];
            
            // track total energies in J = N-m
            //	mp is g, stored energy is N/m^2 cm^3/g, vel is mm/sec
            // workEnergy in J =  1.0e-6*mp*mpm[p]->GetWorkEnergy()
            // plastic 1.0e-6*mp*mpm[p]->GetPlastEnergy()
            // external work 1.e-9*mpm[p]->GetExtWork()
            // kinetic energy 0.5e-9*mp*(vel.x*vel.x+vel.y*vel.y)
            
            // plastic energy per unit thickness (units of N) (only needed energy balance crack growth)
            double mp = mpnt->mp;
            totalPlastic += 1.0e-3*mp*mpnt->GetPlastEnergy()/mpnt->thickness();
            //totalPotential += 1.0e-3*(mp*mpnt->GetStrainEnergy()
            //                        + 0.5e-3*mp*(mpnt->vel.x*mpnt->vel.x+mpnt->vel.y*mpnt->vel.y)
            //                        - 1.e-3*mpnt->GetExtWork())/mpnt->thickness();
			throw "external work is no longer available";
        }
    }
    
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
            {	// crack tip terms last two times propagated
                crkTip->potential[2]=crkTip->potential[1];
                crkTip->potential[1]=crkTip->potential[0];
                crkTip->potential[0]=totalPotential;
                crkTip->plastic[2]=crkTip->plastic[1];
                crkTip->plastic[1]=crkTip->plastic[0];
                crkTip->plastic[0]=totalPlastic;
                crkTip->clength[2]=crkTip->clength[1];
                crkTip->clength[1]=crkTip->clength[0];
                crkTip->clength[0]=nextCrack->Length();
            
                // see if it grows
				int shouldGo=theMaterials[inMat-1]->ShouldPropagate(crkTip,tipDir,nextCrack,fmobj->np,0);
				isAlt[0] = 0;
				if(shouldGo==GROWNOW)
				{	nextCrack->SetAllowAlternate(i,FALSE);
				}
				else if(nextCrack->GetAllowAlternate(i))
				{	shouldGo=theMaterials[inMat-1]->ShouldPropagate(crkTip,tipDir,nextCrack,fmobj->np,1);
					if(shouldGo==GROWNOW) strcpy(isAlt," (alt)");
				}
				if(shouldGo==GROWNOW)
                {   theResult=GROWNOW;
                    tipElem=crkTip->planeInElem-1;
                    if(fabs(tipDir.x)>fabs(tipDir.y))
                        cSize=theElements[tipElem]->xmax-theElements[tipElem]->xmin;
                    else
                        cSize=theElements[tipElem]->ymax-theElements[tipElem]->ymin;
                    grow.x=cellsPerPropagationStep*cSize*tipDir.x;
                    grow.y=cellsPerPropagationStep*cSize*tipDir.y;
					
					// adjust if crossing another crack
					double p = nextCrack->AdjustGrowForCrossing(&grow,crkTip);
					
					// crack number and tip
					archiver->IncrementPropagationCounter();
                    cout << "# propagation" << isAlt << " crack-tip " << nextCrack->GetNumber() << "-" << i;
					
					// summarize
					cout << " at t=" << 1000*mtime << " with J=Jtip+Jzone : " << 1000.*crkTip->Jint.z <<
							" = " << 1000.*crkTip->Jint.x << " + " << 1000.*(crkTip->Jint.z-crkTip->Jint.x) << endl;
                    
                    // if jump is .7 or more cells, make more than 1 segment
                    int iseg,numSegs = 1;
                    if(p*cellsPerPropagationStep>.7) numSegs= 2*(p*cellsPerPropagationStep+.25);
                    CrackSegment *newCrkTip;
                    for(iseg=1;iseg<=numSegs;iseg++)
                    {   growTo.x=crkTip->x+(double)iseg*grow.x/(double)numSegs;
                        growTo.y=crkTip->y+(double)iseg*grow.y/(double)numSegs;
                        if(fmobj->dflag[0]==4) growTo.y=0.;			// force cutting simulation to stay in cut plane at 0
                        newCrkTip=nextCrack->Propagate(growTo,(int)i,theMaterials[inMat-1]->tractionMat[0]);
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
                else
                {   // when no growth, restore previous growth results
                    crkTip->potential[0]=crkTip->potential[1];
                    crkTip->potential[1]=crkTip->potential[2];
                    crkTip->plastic[0]=crkTip->plastic[1];
                    crkTip->plastic[1]=crkTip->plastic[2];
                    crkTip->clength[0]=crkTip->clength[1];
                    crkTip->clength[1]=crkTip->clength[2];
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
void PropagateTask::ArrestGrowth(int newArrest) { arrested=newArrest; }
bool PropagateTask::Arrested(void) { return arrested; }

