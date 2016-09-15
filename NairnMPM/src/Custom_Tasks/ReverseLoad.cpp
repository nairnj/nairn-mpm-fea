/********************************************************************************
    ReverseLoad.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Aug 20 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
	
	To stop or reverse load based on crack length or on critical
        value of a global quantity
	
	Parameters
		crackNum: 0 (default) for any crack on # for crack number (int)
		maxLength: length to change (default 10) (double), for global is is value
		quantity: >=0 to be global quantity instead of crack number
		subcode: used weh quantity is a history variable
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

#pragma mark INITIALIZE

// Constructors
ReverseLoad::ReverseLoad()
{
    crackNum=0;			// means any crack
    finalLength=10.;	// final length in mm
    reversed=CHECKING_PHASE;
    style=REVERSE;
	quantity = -1;
	whichMat = 0;
    holdTime = -1.;
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
    
    // style (out or range will default to REVERSE)
    else if(strcmp(pName,"style")==0)
    {	input=INT_NUM;
        return (char *)&style;
    }
    
	// which material (if !=0) for global quantity
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
			quantity = GlobalQuantity::DecodeGlobalQuantity(quant,&subcode);
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
			sprintf(quant,"%s mat %d",quant,whichMat);
		
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
		cout << quant << " passes " << finalLength << endl;
		
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
		if(propagateTask==NULL)
		{	reversed=REVERSED_PHASE;
			return nextTask;
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
		{	cout << "any crack length reaches " << finalLength
				<< " " << UnitsController::Label(CULENGTH_UNITS) << endl;
		}
		else
		{	cout << "crack " << crackNum << " reaches "
				<< finalLength << " " << UnitsController::Label(CULENGTH_UNITS) << endl;
		}
	}
    
    return nextTask;
}

// Called when custom tasks are all done on a step
// throws CommonException() to end simulation
CustomTask *ReverseLoad::FinishForStep(void)
{
    // change to true trigger or on end of holding phase
	bool status = false;
    
    // after checking decide if stop or if final hold phase
    switch(reversed)
    {   case REVERSED_PHASE:
            // once reversed, check if reached final time and then exit
            // othersize just continue
            if(style==REVERSE && mtime>finalTime)
                throw CommonException("Load or displacement has returned to zero","ReverseLoad::FinishForStep");
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
                status = archiver->PassedLastArchived(quantity,finalLength);
            }
            
            // check cracks
            else if(propagateTask->theResult!=NOGROWTH)
            {	// check desired crack or all cracks
                int cnum=0;
                CrackHeader *nextCrack=firstCrack;
                while(nextCrack!=NULL)
                {	cnum++;
                    if(crackNum==0 || crackNum==cnum)
                    {   // found crack, check the length
                        if(nextCrack->Length()>finalLength) status = true;
                        
                        // exit if critical or only one crack to check
                        if(status || crackNum==cnum) break;
                    }
                    nextCrack=(CrackHeader *)nextCrack->GetNextObject();
                }
            }
            break;
    }
    
    // if status is false, then nothing to do
    if(!status) return nextTask;
	
    // Print message that task has been triggered
    if(reversed == CHECKING_PHASE)
    {   if(quantity>=0)
        {	if(style==ABORT && holdTime<0.)
                throw CommonException("Global quantity has reached specified value","ReverseLoad::FinishForStep");
            
            cout << "# Critical global quantity reached at time t: " << mtime*UnitsController::Scaling(1.e3)
				<< " " << UnitsController::Label(ALTTIME_UNITS) << endl;
        }
        else
        {	if(style==ABORT && holdTime<0.)
                throw CommonException("Crack has reached specified length","ReverseLoad::FinishForStep");
            
            // stop propgation
            propagateTask->ArrestGrowth(TRUE);
            cout << "# Crack growth arrested at time t: " << mtime*UnitsController::Scaling(1.e3)
				<< " " << UnitsController::Label(ALTTIME_UNITS) << endl;
        }
        
        // final time if reversing
        finalTime = 2.*mtime;
        
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
            throw CommonException("ReverseLoad task hold time has ended","ReverseLoad::FinishForStep");
        
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
        
        // reverse rigid particles at constant velocity
        int p;
        for(p=nmpmsNR;p<nmpms;p++)
        {	RigidMaterial *mat = (RigidMaterial *)theMaterials[mpm[p]->MatID()];
            if(mat->IsConstantVelocity())
            {	if(style==REVERSE)
                    mpm[p]->ReverseParticle(reversed==HOLDING_PHASE,holdTime>0.);
                else
                    mpm[p]->StopParticle();
            }
        }
    }
	
	return nextTask;
}

