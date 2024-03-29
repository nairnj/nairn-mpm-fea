/********************************************************************************
    ReverseLoad.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Aug 20 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
	
	To stop or reverse load based on crack length or on critical
        value of a global quantity
	
	Parameters
		crackNum: 0 (default) for any crack or # for crack number (int)
		maxLength: length to change (default 10) (double), for global is is value
		quantity: >=0 to be global quantity instead of crack number
		subcode: used when quantity is a history variable or crack number
		whichMat: used when quantity is for one material type
		maxValue: max value (stored in maxLength)
		style: What to do when maxLength is reached
			0: REVERSE (default): Reverse all linear loads and stop for good when
					loads drop to zero again. Also reverse all rigid particles.
			1: HOLD; Tell all linear loads to stay constant and continue analysis.
					Also stop all rigid particles
			2: NOCHANGE: Crack stopped, but loading continued
			3: ABORT: Stop the analysis immediately
********************************************************************************/

#include "stdafx.h"
#include "Custom_Tasks/ReverseLoad.hpp"
#include "Custom_Tasks/PropagateTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Boundary_Conditions/MatPtLoadBC.hpp"
#include "Boundary_Conditions/MatPtTractionBC.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/RigidMaterial.hpp"
#include "Exceptions/CommonException.hpp"
#include "Global_Quantities/GlobalQuantity.hpp"
#include "System/ArchiveData.hpp"
#include "System/UnitsController.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"

#pragma mark INITIALIZE

// Constructors
ReverseLoad::ReverseLoad() : CustomTask()
{
    crackNum=0;			// means any crack
    finalLength=10.;	// final length in mm
    reversed=CHECKING_PHASE;
    style=REVERSE;
	quantity = -1;
	whichMat = 0;
    holdTime = -1.;
    debondLength = 0;   // made true if no propagation tasks
    waitForDrop = false;
    maxValue = 0.;
    passedMinValue = false;
}

// Return name of this task
const char *ReverseLoad::TaskName(void) { return "Reverse Load Task"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
char *ReverseLoad::InputParam(char *pName,int &input,double &gScaling)
{
    // crack number to check or 0 to check all cracks
    if(strcmp(pName,"crackNumber")==0)
    {	input=INT_NUM;
        return (char *)&crackNum;
    }
    
    // crack length (in mm) or critical value if task triggered on global quantity
    else if(strcmp(pName,"maxLength")==0 || strcmp(pName,"maxValue")==0)
    {	input=DOUBLE_NUM;
        return (char *)&finalLength;
    }
    
    // minimum value for peak monitoring
    else if(strcmp(pName,"minValue")==0)
    {   input=DOUBLE_NUM;
        waitForDrop = true;
        return (char *)&minValue;
    }
    
   if(strcmp(pName,"debondedLength")==0)
    {   input=INT_NUM;
        return (char *)&debondLength;
    }
    
    // style (out of range will default to REVERSE)
    else if(strcmp(pName,"style")==0)
    {	input=INT_NUM;
        return (char *)&style;
    }
    
	// which material (if !=0), but only for global quantity
    else if(strcmp(pName,"mat")==0)
    {	input=INT_NUM;
        return (char *)&whichMat;
    }
	
	// hold time
    else if(strcmp(pName,"hold")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&holdTime,gScaling,1.e-3);
    }
	
	// look for "global quantity" instead of crack growth
    else if(strlen(pName)>7)
	{	strcpy(quant,pName);
		quant[6]=0;
		if(strcmp(quant,"global")==0)
		{	strcpy(quant,&pName[7]);
			quantity = GlobalQuantity::DecodeGlobalQuantity(quant,&subcode,&whichMat);
			input=INT_NUM;
			return (char *)&ignoreArgument;
		}
    }
	
	return CustomTask::InputParam(pName,input,gScaling);
}

#pragma mark GENERIC TASK METHODS

// at beginning of analysis
// throws CommonException()
CustomTask *ReverseLoad::Initialize(void)
{
	// is it for a global quantity?
	if(quantity>=0)
	{	// finish name
		if(whichMat!=0)
			snprintf(quant,quantSize,"%s mat %d",quant,whichMat);
		
		// find the global quantity
		switch(style)
		{	case HOLD:
				cout << "Loading held if ";
                if(holdTime>0.) holdTime = -1.;
				break;
			case NOCHANGE:
				throw CommonException("A ReverseLoad task cannot continue with no change when based on global archive value","ReverseLoad::Initialize");
				break;
			case ABORT:
                if(holdTime>0.)
				{	cout << "Analysis held " << holdTime*UnitsController::Scaling(1.e3)
						<< " " << UnitsController::Label(ALTTIME_UNITS) << " and then stopped if ";
				}
                else
                    cout << "Analysis stopped if ";
				break;
			default:
				style=REVERSE;			// make all others equal to REVERSE
                if(holdTime>0.)
				{	cout << "Analysis held " << holdTime*UnitsController::Scaling(1.e3)
						<< " " << UnitsController::Label(ALTTIME_UNITS) << " and then loading reversed to zero if ";
				}
                else
                    cout << "Loading reversed to zero if ";
				break;
		}
        if(!waitForDrop)
            cout << quant << " passes " << finalLength << endl;
        else
        {   cout << quant << " passes " << minValue << endl;
            cout << "   and then drops to " << finalLength << " of peak value" << endl;
            if(finalLength<=0. || finalLength>=1.)
            {   throw CommonException("For peak tracking, the maxValue must be between 0 and 1","ReverseLoad::Initialize");
            }
       }
		
		// convert to index
		int qIndex=0;
		GlobalQuantity *nextGlobal=firstGlobal;
		while(nextGlobal!=NULL)
		{	if(nextGlobal->IsSameQuantity(quantity,subcode,whichMat)) break;
			qIndex++;
			nextGlobal=nextGlobal->GetNextGlobal();
		}
		
		// error if not found
		if(nextGlobal==NULL)
			throw CommonException("A ReverseLoad task's global quantity is not available in the current global archives","ReverseLoad::Initialize");
		
		// store in quantity
		quantity = qIndex;
	}
	
	else
    {	// skipped if no propagation being done
        // No runs even if no propgation, but only debond growth makes sense
		if(propagateTask==NULL)
		{	//reversed=REVERSED_PHASE;
			//return nextTask;
            debondLength = 1;
		}
		
		switch(style)
		{	case HOLD:
				cout << "Crack stopped and load held if ";
                if(holdTime>0.) holdTime = -1.;
				break;
			case NOCHANGE:
				cout << "Crack stopped and load continued if ";
                if(holdTime>0.) holdTime = -1.;
				break;
			case ABORT:
                if(holdTime>0.)
				{	cout << "Cracks stopped and analysis stopped after hold of " << holdTime*UnitsController::Scaling(1.e3)
						<< " " << UnitsController::Label(ALTTIME_UNITS) << " if ";
				}
                else
                    cout << "Cracks and analysis stopped if ";
				break;
			default:
				style=REVERSE;			// make all others equal to REVERSE
                if(holdTime>0.)
				{	cout << "Cracks stopped and loads reversed to zero after hold of " << holdTime*UnitsController::Scaling(1.e3)
						<< " " << UnitsController::Label(ALTTIME_UNITS) << " if ";
				}
                else
                    cout << "Cracks stopped and loads reversed to zero if ";
				break;
		}
		
		if(crackNum==0)
        {	cout << "any crack's";
            if(debondLength!=0) cout << " debonded";
            cout <<" length reaches " << finalLength << " " << UnitsController::Label(CULENGTH_UNITS) << endl;
		}
		else
        {	cout << "crack " << crackNum << "'s";
            if(debondLength!=0) cout << " debonded";
            cout << " length reaches " << finalLength << " " << UnitsController::Label(CULENGTH_UNITS) << endl;
		}
	}
    
    return nextTask;
}

// Called when custom tasks are all done on a step
// throws CommonException() to end simulation
CustomTask *ReverseLoad::FinishForStep(bool &removeMe)
{
    // change to true if triggered or on end of holding phase
	bool status = false;
    
    // after checking decide if stop or if final hold phase
    switch(reversed)
    {   case REVERSED_PHASE:
            // once reversed, check if reached final time and then exit
            // othersize just continue
            if(style==REVERSE && mtime>finalTime)
                throw CommonException("Load or displacement has returned to zero","ReverseLoad::FinishForStep",noErr);
            break;
        
        case HOLDING_PHASE:
            // change status to true if done with holding phase
            if(mtime>endHoldTime) status=true;
            break;
        
        case CHECKING_PHASE:
        default:
            // check global quantity
            if(quantity>=0)
            {	// see if has passed desider value
                double currentValue = archiver->GetLastArchived(quantity);
                
                if(!waitForDrop)
                {   // just check if passed specified value
                    if(finalLength>=0.)
                        status = (currentValue >= finalLength);
                    else
                        status = (currentValue <= finalLength);
                }
                else if(minValue>=0.)
                {   // looking for drop to finalLength*maxValue after passing minValue
                    if(currentValue>maxValue) maxValue = currentValue;
                    if(!passedMinValue && maxValue>=minValue) passedMinValue = true;
                    status = false;
                    if(passedMinValue)
                    {   if(currentValue<0.)
                        {   // went nagative
                            status = true;
                        }
                        else if(currentValue/maxValue < finalLength)
                            status = true;
                    }
                 }
                else
                {   // looking for negative results (pass -minValue then drop to fraction of maxValue)
                    if(currentValue<maxValue) maxValue = currentValue;
                    if(!passedMinValue && maxValue<=minValue) passedMinValue = true;
                    status = false;
                    if(passedMinValue)
                    {   if(currentValue>0.)
                        {   // went positive
                            status = true;
                        }
                        else if(currentValue/maxValue < finalLength)
                            status = true;
                    }
                 }
            }

            // check cracks
            else
            {   bool checkThisStep = propagateTask!=NULL ? propagateTask->theResult!=NOGROWTH : true ;
                if(debondLength!=0) checkThisStep = true;         // always check if might just be debonds
                if(checkThisStep)
                {	// check desired crack or all cracks
                    int cnum=0;
                    CrackHeader *nextCrack=firstCrack;
                    while(nextCrack!=NULL)
                    {	cnum++;
                        if(crackNum==0 || crackNum==cnum)
                        {   // found crack, check the length
                            double crackLength = debondLength!=0 ? nextCrack->DebondedLength() : nextCrack->Length() ;
                            if(crackLength>finalLength) status = true;
                            
                            // exit if critical or only one crack to check
                            if(status || crackNum==cnum) break;
                        }
                        nextCrack=(CrackHeader *)nextCrack->GetNextObject();
                    }
                }
            }
            break;
    }
    
    // if status is false, then nothing to do
    if(!status) return nextTask;
	
    // Print message that task has been triggered
	// throw is aborting, or message line if continuing
    if(reversed == CHECKING_PHASE)
    {   if(quantity>=0)
        {	if(style==ABORT && holdTime<0.)
			{	char rev[200];
                size_t revSize=200;
				if(waitForDrop)
					snprintf(rev,revSize,"Value %s reached %.6g and then dropped to %.6g",quant,maxValue,finalLength*maxValue);
				else
                    snprintf(rev,revSize,"Value %s reached %.6g",quant,finalLength);
				throw CommonException(rev,"ReverseLoad::FinishForStep",noErr);
			}
			
			if(waitForDrop)
			{	cout << "# " << quant << " reached " << maxValue
						<< " and then dropped to " << (finalLength*maxValue);
			}
			else
			{	cout << "# " << quant << " reached " << finalLength;
			}
			cout << " at time t: " << mtime*UnitsController::Scaling(1.e3)
				<< " " << UnitsController::Label(ALTTIME_UNITS) << endl;
        }
        else
        {	if(style==ABORT && holdTime<0.)
			{	char rev[200];
                size_t revSize=200;
				snprintf(rev,revSize,"Crack has reached length = %.6g",finalLength);
                throw CommonException(rev,"ReverseLoad::FinishForStep",noErr);
			}
            
            // stop propgation
            if(propagateTask!=NULL)
            {   propagateTask->ArrestGrowth(true);
                cout << "# Crack length = " << finalLength << " and growth arrested at time t: "
							<< mtime*UnitsController::Scaling(1.e3)
                            << " " << UnitsController::Label(ALTTIME_UNITS) << endl;
            }
        }
        
        // final time if reversing
        
        // if holding, wait to trigger change until latter
        if(holdTime>0.)
        {   reversed = HOLDING_PHASE;
            endHoldTime = mtime + holdTime;
            finalTime = 2.*mtime + holdTime;
        }
        else
        {   reversed = REVERSED_PHASE ;
            finalTime = 2.*mtime;
        }
    }
    else
    {   // if an abort after holding, exit now
       if(style==ABORT)
            throw CommonException("ReverseLoad task hold time has ended","ReverseLoad::FinishForStep",noErr);
        
        // message that hold time is over
        cout << "# ReverseLoad task hold time ended at t:" << mtime*UnitsController::Scaling(1.e3)
			<< " " << UnitsController::Label(ALTTIME_UNITS) << endl;
        reversed = REVERSED_PHASE;
    }
    
    // REVERSE: reverse linear loads, zero constant loads, reverse or zero constant velocity rigid particles
    // HOLD: make linear loads constant, stop constant velocity rigid particles
    // ABORT: here is holding, treat like a HOLD
    // finalTime will be twice current time, or if any reversed linear loads, the time last one gets to zero
    if(style!=NOCHANGE)
    {   // load BCs
        MatPtLoadBC *nextLoad=firstLoadedPt;
        while(nextLoad!=NULL)
        {	if(style==REVERSE)
                nextLoad=nextLoad->ReverseLinearLoad(mtime,&finalTime,reversed==HOLDING_PHASE);
            else
                nextLoad=nextLoad->MakeConstantLoad(mtime);
        }
        
        // traction BCs
        MatPtTractionBC *nextTraction=firstTractionPt;
        while(nextTraction!=NULL)
        {	if(style==REVERSE)
                nextTraction=(MatPtTractionBC *)nextTraction->ReverseLinearLoad(mtime,&finalTime,reversed==HOLDING_PHASE);
            else
                nextTraction=(MatPtTractionBC *)nextTraction->MakeConstantLoad(mtime);
        }
        
        // reverse rigid contact and BC particles
        int p;
        for(p=nmpmsRB;p<nmpms;p++)
		{	if(mpm[p]->InReservoir()) continue;
			RigidMaterial *mat = (RigidMaterial *)theMaterials[mpm[p]->MatID()];
            if(mat->IsConstantVelocity())
            {	if(style==REVERSE)
                    mpm[p]->ReverseParticle(reversed==HOLDING_PHASE,holdTime>0.);
                else
                    mpm[p]->StopParticle();
            }
        }
		
		// velocity BCs set to zero, will not reverse
		NodalVelBC::holdAllVelocityBCs = true;
    }
	
	return nextTask;
}

