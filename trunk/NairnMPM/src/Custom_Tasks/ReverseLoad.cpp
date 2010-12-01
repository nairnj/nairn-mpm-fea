/********************************************************************************
    ReverseLoad.cpp
    NairnMPM
    
    Created by John Nairn on Wed Aug 20 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
	
	To stop or reverse load based on crack length
	
	Parameters
		crackNum: 0 (default) for any crack on # for crack number (int)
		maxLength: length to change (default 10) (double)
		style: What to do when maxLength is reached
			0: REVERSE (default): Reverse all linear loads and stop for good when
					loads drop to zero again. Also reverse all rigid particles.
			1: HOLD; Tell all linear loads to stay constant and continue analysis.
					Also stop all rigid particles
			2: NOCHANGE: Crack stopped, but loading continued
			3: ABORT: Stop the analysis immediately
********************************************************************************/

#include "Custom_Tasks/ReverseLoad.hpp"
#include "Custom_Tasks/PropagateTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Boundary_Conditions/MatPtLoadBC.hpp"
#include "Exceptions/MPMTermination.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"

#pragma mark INITIALIZE

// Constructors
ReverseLoad::ReverseLoad()
{
    crackNum=0;			// means any crack
    finalLength=10.;	// final length in mm
    reversed=FALSE;
    style=REVERSE;
}

// Return name of this task
const char *ReverseLoad::TaskName(void) { return "Reverse Load Task"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
char *ReverseLoad::InputParam(char *pName,int &input)
{
    // crack number to check or 0 to check all cracks
    if(strcmp(pName,"crackNumber")==0)
    {	input=INT_NUM;
        return (char *)&crackNum;
    }
    
    // crack length in mm
    else if(strcmp(pName,"maxLength")==0)
    {	input=DOUBLE_NUM;
        return (char *)&finalLength;
    }
    
    // style (out or range will default to REVERSE)
    else if(strcmp(pName,"style")==0)
    {	input=INT_NUM;
        return (char *)&style;
    }
    
    return CustomTask::InputParam(pName,input);
}

#pragma mark GENERIC TASK METHODS

// at beginning of analysis
CustomTask *ReverseLoad::Initialize(void)
{
    // skipped if no propagation being done
    if(propagateTask==NULL)
    {	reversed=TRUE;
        return nextTask;
    }
    
    switch(style)
    {	case HOLD:
            cout << "Crack stopped and load held if ";
            break;
        case NOCHANGE:
            cout << "Crack stopped and load continued if ";
            break;
		case ABORT:
			cout << "Crack and analysis stopped if ";
            break;
		default:
            style=REVERSE;			// make all others equal to REVERSE
            cout << "Crack stopped and load reversed to zero if ";
            break;
    }
    
    if(crackNum==0)
    {	cout << "any crack length reaches " << finalLength
            << " mm" << endl;
    }
    else
    {	cout << "crack " << crackNum << " reaches "
            << finalLength << " mm" << endl;
    }
    return nextTask;
}

// Called when custom tasks are all done on a step
CustomTask *ReverseLoad::FinishForStep(void)
{
    CrackHeader *nextCrack;
    int cnum;
    MatPtLoadBC *nextLoad;
    double currentLoad,loadx,loady;
    
    // exit if done or no need
    if(reversed)
	{	if(style==REVERSE)
		{	if(firstLoadedPt!=NULL)
			{	Vector *pFext=mpm[(firstLoadedPt->ptNum)-1]->GetPFext();
				loadx=pFext->x;
				loady=pFext->y;
				currentLoad=fmax(fabs(loadx),fabs(loady));
				if(loadx<0. || loady<0.) currentLoad=-currentLoad;
				if(currentLoad*finalLoad<0.)
					throw MPMTermination("Load has returned to zero","ReverseLoad::FinishForStep");
			}
			else if(mtime>=finalLoad)
				throw MPMTermination("Displacement has returned to zero","ReverseLoad::FinishForStep");
		}
		return nextTask;
	}
	
	// exit if not propagating
    if(propagateTask->theResult==NOGROWTH) return nextTask;
    
    // check desired crack or all cracks
    cnum=0;
    nextCrack=firstCrack;
    while(nextCrack!=NULL)
    {	cnum++;
        if(crackNum==0 || crackNum==cnum)
        {   // check the length
            if(nextCrack->Length()>finalLength)
			{	// is ABORT, then exit analysis now
				if(style==ABORT)
					throw MPMTermination("Crack has reached specified length","ReverseLoad::FinishForStep");
					
            	// stop propgation
                propagateTask->ArrestGrowth(TRUE);
                cout << "# Crack arrested at time t: " << 1000* mtime << endl;

                // reverse the loads and rigid particles (REVERSE or HOLD)
                if(style!=NOCHANGE)
                {   reversed=TRUE;
                    nextLoad=firstLoadedPt;
                    while(nextLoad!=NULL)
                    {	if(style==REVERSE)
                            nextLoad=nextLoad->ReverseLinearLoad(mtime);
                        else
                            nextLoad=nextLoad->MakeConstantLoad(mtime);
                    }
					
					// reverse rigid particles
					int p;
					for(p=0;p<nmpms;p++)
					{	if(theMaterials[mpm[p]->MatID()]->Rigid())
						{	if(style==REVERSE)
								mpm[p]->ReverseParticle();
							else
								mpm[p]->StopParticle();
						}
					}
					
                }
				
				// find when to stop based on first point
				if(style==REVERSE)
				{	if(firstLoadedPt!=NULL)
					{	Vector *pFext=mpm[(firstLoadedPt->ptNum)-1]->GetPFext();
						loadx=pFext->x;
						loady=pFext->y;
						finalLoad=fmax(fabs(loadx),fabs(loady));
						if(loadx<0. || loady<0.) finalLoad=-finalLoad;
					}
					else
						finalLoad=2.*mtime;
				}
            }
        }
        nextCrack=(CrackHeader *)nextCrack->GetNextObject();
    }
    
    return nextTask;
}

