/********************************************************************************
    ReverseLoad.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Aug 20 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
	
	To stop or reverse load based on crack length
	
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

#pragma mark INITIALIZE

int ignoreArg;

// Constructors
ReverseLoad::ReverseLoad()
{
    crackNum=0;			// means any crack
    finalLength=10.;	// final length in mm
    reversed=FALSE;
    style=REVERSE;
	quantity = -1;
	whichMat = 0;
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
    else if(strcmp(pName,"maxLength")==0 || strcmp(pName,"maxValue")==0)
    {	input=DOUBLE_NUM;
        return (char *)&finalLength;
    }
    
    // style (out or range will default to REVERSE)
    else if(strcmp(pName,"style")==0)
    {	input=INT_NUM;
        return (char *)&style;
    }
    
	// convert to reverse on reaching critical value of a quantity
    else if(strcmp(pName,"mat")==0)
    {	input=INT_NUM;
        return (char *)&whichMat;
    }
	
	// look for "global quantity"
    else if(strlen(pName)>7)
	{	strcpy(quant,pName);
		quant[6]=0;
		if(strcmp(quant,"global")==0)
		{	strcpy(quant,&pName[7]);
			quantity = GlobalQuantity::DecodeGlobalQuantity(quant,&subcode);
			input=INT_NUM;
			return (char *)&ignoreArg;
		}
    }
	
	return CustomTask::InputParam(pName,input);
}

#pragma mark GENERIC TASK METHODS

// at beginning of analysis
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
				break;
			case NOCHANGE:
				throw CommonException("The ReverseLoad cannot continue with no change when based on global archive value","ReverseLoad::Initialize");
				break;
			case ABORT:
				cout << "Analysis stopped if ";
				break;
			default:
				style=REVERSE;			// make all others equal to REVERSE
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
			throw CommonException("The ReverseLoad global quantity is not available in the current global archives","ReverseLoad::Initialize");
		
		// store in quantity
		quantity = qIndex;
	}
	
	else
    {	// skipped if no propagation being done
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
	}
	
    return nextTask;
}

// Called when custom tasks are all done on a step
// throw CommonException() to end simulation
CustomTask *ReverseLoad::FinishForStep(void)
{
    CrackHeader *nextCrack;
    int cnum;
    MatPtLoadBC *nextLoad;
	MatPtTractionBC *nextTraction;
    
    // exit if done or no need
    if(reversed)
	{	if(style==REVERSE && mtime>finalTime)
			throw CommonException("Load or displacement has returned to zero","ReverseLoad::FinishForStep");
		return nextTask;
	}
	
	// change to true if critical
	bool status = false;
	
	// check global quantity
	if(quantity>=0)
	{	// see if has passed desider value
		status = archiver->PassedLastArchived(quantity,finalLength);
	}
	
	// check cracks
	else if(propagateTask->theResult!=NOGROWTH)
	{	// check desired crack or all cracks
		cnum=0;
		nextCrack=firstCrack;
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
    
	// Reverse load if needed
	if(status)
	{	// is ABORT, then exit analysis now
		if(quantity>=0)
		{	if(style==ABORT)
				throw CommonException("Global quantity has reached specified value","ReverseLoad::FinishForStep");
			
			cout << "# Critical global quantity reached at time t: " << 1000* mtime << endl;
		}
		else
		{	if(style==ABORT)
				throw CommonException("Crack has reached specified length","ReverseLoad::FinishForStep");
			
			// stop propgation
			propagateTask->ArrestGrowth(TRUE);
			cout << "# Crack arrested at time t: " << 1000* mtime << endl;
		}
		
		// REVERSE: reverse linear loads, zero constant loads, reverse or zero constant velocity rigid particles
		// HOLD: make linear loads constant, stop constant velocity rigid particles
		// finalTime will be twice current time, or if any reversed linear loads, the time last one gets to zero
		finalTime = 2.*mtime;
		if(style!=NOCHANGE)
		{   reversed=TRUE;
			
			// load BCs
			nextLoad=firstLoadedPt;
			while(nextLoad!=NULL)
			{	if(style==REVERSE)
					nextLoad=nextLoad->ReverseLinearLoad(mtime,&finalTime);
				else
					nextLoad=nextLoad->MakeConstantLoad(mtime);
			}
			
			// traction BCs
			nextTraction=firstTractionPt;
			while(nextTraction!=NULL)
			{	if(style==REVERSE)
					nextTraction=(MatPtTractionBC *)nextTraction->ReverseLinearLoad(mtime,&finalTime);
				else
					nextTraction=(MatPtTractionBC *)nextTraction->MakeConstantLoad(mtime);
			}
			
			// reverse rigid particles at constant velocity
			int p;
			for(p=nmpmsNR;p<nmpms;p++)
			{	RigidMaterial *mat = (RigidMaterial *)theMaterials[mpm[p]->MatID()];
				if(mat->IsConstantVelocity())
				{	if(style==REVERSE)
						mpm[p]->ReverseParticle();
					else
						mpm[p]->StopParticle();
				}
			}
		}
	}
	
	return nextTask;
}

