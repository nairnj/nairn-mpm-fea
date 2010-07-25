/********************************************************************************
    PropagateTask.cpp
    NairnMPM
    
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

// globals
PropagateTask *propagateTask=NULL;
double PropagateTask::cellsPerPropagationStep=0.5;

#pragma mark Constructors and Destructors

// Constructors
PropagateTask::PropagateTask()
{
    propagateTask=this;
    arrested=FALSE;
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
CustomTask *PropagateTask::PrepareForStep(bool &doPropExtraps)
{
    theResult=NOGROWTH;
    
    if(mtime+timestep>=nextPropTime && !arrested)
    {	// ready to do crack propagation
    	doPropCalcs=doPropExtraps=TRUE;
        nextPropTime+=propTime;
        
        // Make sure J and K are available if needed for current criterion
		CrackHeader *nextCrack=firstCrack;
		while(nextCrack!=NULL)
		{	theJKTask->ScheduleJK(nextCrack->CriterionNeeds());
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
    CrackHeader *nextCrack;
    CrackSegment *crkTip;
    double cSize;
    int i,inMat;
    Vector tipDir,growTo,grow;
    int tipElem;
    
    // if not needed, just exit
    if(!doPropCalcs) return nextTask;
    
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
                if(theMaterials[inMat-1]->ShouldPropagate(crkTip,tipDir,nextCrack,fmobj->np)==GROWNOW)
                {   theResult=GROWNOW;
                    tipElem=crkTip->planeInElem-1;
                    if(fabs(tipDir.x)>fabs(tipDir.y))
                        cSize=theElements[tipElem]->xmax-theElements[tipElem]->xmin;
                    else
                        cSize=theElements[tipElem]->ymax-theElements[tipElem]->ymin;
                    grow.x=cellsPerPropagationStep*cSize*tipDir.x;
                    grow.y=cellsPerPropagationStep*cSize*tipDir.y;
                    growTo.x=crkTip->x+grow.x;
                    growTo.y=crkTip->y+grow.y;
					cout << "# propagation at t=" << 1000*mtime << " with J=Jtip+Jzone : " << crkTip->Jint.z <<
							" = " << crkTip->Jint.x << " + " << crkTip->Jint.z-crkTip->Jint.x << endl;
                    crkTip=nextCrack->Propagate(growTo,(int)i,theMaterials[inMat-1]->tractionMat);
                    
					if(crkTip!=NULL)
					{	// check if crack speed is being controlled
						if(theMaterials[inMat-1]->ControlCrackSpeed(crkTip,propTime))
							nextPropTime=mtime+propTime;
						
						// crack tip heating (if activated)
						if(ConductionTask::active)
							conduction->StartCrackTipHeating(crkTip,grow,nextCrack->thickness);
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

#pragma mark TASK EXTRAPOLATION METHODS

// initialize for crack extrapolations
CustomTask *PropagateTask::BeginExtrapolations(void)
{
    // skip if already set up
    if(!doPropCalcs) return nextTask;

    totalPlastic=0.;
    totalPotential=0.;
    
    return nextTask;
}

// add particle data to some calculation
CustomTask *PropagateTask::ParticleExtrapolation(MPMBase *mpnt)
{
    double mp;
    
    // skip if not needed
    if(!doPropCalcs) return nextTask;
    
    // track total energies in J = N-m
    //	mp is g, stored energy is N/m^2 cm^3/g, vel is mm/sec
    // strainEnergy 1.0e-6*mp*mpm[p]->GetStrainEnergy()
    // plastic 1.0e-6*mp*mpm[p]->GetPlastEnergy()
    // external work 1.e-9*mpm[p]->GetExtWork()
    // kinetic energy 0.5e-9*mp*(vel.x*vel.x+vel.y*vel.y)
    
    // plastic energy per unit thickness (units of N) (only needed in some crack growth)
    mp=mpnt->mp;
    totalPlastic+=1.0e-3*mp*mpnt->GetPlastEnergy()/mpnt->thickness();
    totalPotential+=1.0e-3*(mp*mpnt->GetStrainEnergy()
        + 0.5e-3*mp*(mpnt->vel.x*mpnt->vel.x+mpnt->vel.y*mpnt->vel.y)
        - 1.e-3*mpnt->GetExtWork())/mpnt->thickness();
    
    return nextTask;
}


