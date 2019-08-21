/********************************************************************************
	Periodic XPIC Custon Task
	nairn-mpm-fea
 
	Created by John Nairn on 1/26/18.
	Copyright (c) 2018 John A. Nairn, All rights reserved.
 
	need method for main code to know this task is active
 ********************************************************************************/

#include "stdafx.h"
#include "Custom_Tasks/PeriodicXPIC.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "System/UnitsController.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Cracks/CrackHeader.hpp"

#pragma mark Constructors and Destructors

// Constructors
PeriodicXPIC::PeriodicXPIC()
{
	verbose = 0;
	periodicXPICorder = 0;
	
	periodicTime = -1.;
	periodicCFL = -1.;
	periodicSteps = -1;
}

// Return name of this task
const char *PeriodicXPIC::TaskName(void) { return "Periodic XPIC Implementation"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
char *PeriodicXPIC::InputParam(char *pName,int &input,double &gScaling)
{
	if(strcmp(pName,"verbose")==0)
	{	input=INT_NUM;
		return (char *)&verbose;
	}
	
	else if(strcmp(pName,"XPICOrder")==0)
	{	input=INT_NUM;
		return (char *)&periodicXPICorder;
	}
	
	// mechanics
	else if(strcmp(pName,"periodicTime")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&periodicTime,gScaling,1.e-3);
	}

	else if(strcmp(pName,"periodicCFL")==0)
	{	input=DOUBLE_NUM;
		return (char *)&periodicCFL;
	}
	
	else if(strcmp(pName,"periodicSteps")==0)
	{	input=INT_NUM;
		return (char *)&periodicSteps;
	}
	
	
	// check remaining commands
	return CustomTask::InputParam(pName,input,gScaling);
}

// Called while reading when done setting all parameters
// Make and needed settings, throw SAXExecption on error
void PeriodicXPIC::Finalize(void)
{
	if(periodicXPICorder<1)
		ThrowSAXException("PeriodicXPIC custom task must set XPIC order to 1 or higher");
	
	// steps override times and chack for active
	bool xpicActive = false;
	
	// mechanics
	if(periodicSteps>0)
	{	periodicTime = periodicCFL = -1.;
		xpicActive = true;
	}
	else if(periodicTime>0. || periodicCFL>0.)
		xpicActive = true;

	// error if nothing
	if(!xpicActive)
		ThrowSAXException("PeriodicXPIC custom task did not activate any XPIC or FMPM calculations");
	
	if(periodicXPICorder>1)
	{	// this will change Preliminary calcs if needed for contact
		bodyFrc.SetXPICVectors(3);
	}
	else
		bodyFrc.SetXPICVectors(1);
}

// called once at start of MPM analysis - initialize and print info
CustomTask *PeriodicXPIC::Initialize(void)
{
	cout << "Periodically switch between " << GetType() << " and FLIP calculations." << endl;
	cout << "   " << GetType() << " order k = " << periodicXPICorder << endl;

	// Mechanics XPIC
	cout << "   Periodic " << GetType() << " interval for mechanics: ";
	if(periodicTime>0. || periodicCFL>0.)
	{	if(periodicCFL>0.) periodicTime = periodicCFL*fmobj->timeStepMinMechanics;
		cout << periodicTime*UnitsController::Scaling(1.e3) << " " << UnitsController::Label(ALTTIME_UNITS) << endl;
		nextPeriodicTime = periodicTime;
		if(periodicTime<=timestep)
		{	periodicTime = -1.;
			periodicSteps = 1;
		}
	}
	else if(periodicSteps>0)
	{	if(periodicSteps>1)
			cout << periodicSteps << " time steps" << endl;
		else
			cout << "every time step" << endl;
		nextPeriodicStep = periodicSteps;
	}
	else
		cout << "not used" << endl;
	
	// when used
	if(periodicTime>0. || periodicSteps>0)
	{	// velocity only not allowed in XPIC
#if MM_XPIC == 1
		if(bodyFrc.XPICVectors()>3)
			cout << "      " << GetType() << " accounts for contact" << endl;
#endif
	}

	// Run for first time step
	bool needExtraps;
	PrepareForStep(needExtraps);
	StepCalculation();
	
	return nextTask;
}

#pragma mark GENERIC TASK METHODS

// called when MPM step is getting ready to do custom tasks
// never uses extrapolations so no need to set needExtraps
CustomTask *PeriodicXPIC::PrepareForStep(bool &needExtraps)
{
	doXPIC = false;

#ifndef TRANSPORT_ONLY
	if(periodicTime>0.)
	{	// controlledby step time
		if(mtime+timestep>=nextPeriodicTime)
		{	doXPIC = true;
			nextPeriodicTime += periodicTime;
		}
	}
	else if(periodicSteps>0)
	{	// controlled by step count
		if(fmobj->mstep+1>=nextPeriodicStep)
		{	doXPIC = true;
			nextPeriodicStep += periodicSteps;
		}
	}
#endif

	return nextTask;
}

// Adjust time step now
CustomTask *PeriodicXPIC::StepCalculation(void)
{
	// reset to FLIP
	bodyFrc.SetXPICOrder(0);
	bodyFrc.SetUsingVstar(VSTAR_NOT_USED);
	
#ifndef TRANSPORT_ONLY
	// mechanics XPIC
	if(doXPIC)
	{	// switch to XPIC for next time step
		bodyFrc.SetXPICOrder(periodicXPICorder);
		if(periodicXPICorder>1)
		{	bodyFrc.SetUsingVstar(VSTAR_NO_CONTACT);
#if MM_XPIC == 1
			if(firstCrack!=NULL || fmobj->multiMaterialMode)
				bodyFrc.SetUsingVstar(VSTAR_WITH_CONTACT);
#endif
		}
		else
			bodyFrc.SetUsingVstar(VSTAR_NOT_USED);
	
		if(verbose)
		{	cout <<"# Use " << GetType() << "(" << periodicXPICorder << ") in step " << (fmobj->mstep+1) << endl;
		}
	}
#endif
	
	return nextTask;
}

// return the type
const char *PeriodicXPIC::GetType(void)
{
	return "XPIC";
}
