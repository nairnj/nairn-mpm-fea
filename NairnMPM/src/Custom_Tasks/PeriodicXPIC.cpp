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

	periodicFMPMorder = 0;
	gridBCOption = -1;
	
	periodicTimeConduction = periodicTimeDiffusion = -1.;
	periodicCFLConduction = periodicCFLDiffusion = -1.;
	periodicStepsConduction = periodicStepsDiffusion = -1;
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
	
	// conduction
	else if(strcmp(pName,"periodicTimeConduction")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&periodicTimeConduction,gScaling,1.e-3);
	}
	
	else if(strcmp(pName,"periodicCFLConduction")==0)
	{	input=DOUBLE_NUM;
		return (char *)&periodicCFLConduction;
	}
	
	else if(strcmp(pName,"periodicStepsConduction")==0)
	{	input=INT_NUM;
		return (char *)&periodicStepsConduction;
	}
	
	// diffusion
	else if(strcmp(pName,"periodicTimeDiffusion")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&periodicTimeDiffusion,gScaling,1.e-3);
	}
	
	else if(strcmp(pName,"periodicCFLDiffusion")==0)
	{	input=DOUBLE_NUM;
		return (char *)&periodicCFLDiffusion;
	}
	
	else if(strcmp(pName,"periodicStepsDiffusion")==0)
	{	input=INT_NUM;
		return (char *)&periodicStepsDiffusion;
	}
	
	else if(strcmp(pName,"FMPMOrder")==0)
	{	input=INT_NUM;
		return (char *)&periodicFMPMorder;
	}
	
	else if(strcmp(pName,"GridBCOption")==0)
	{	input=TEXT_PARAMETER;
		return (char *)&gridBCOption;
	}
	
	// check remaining commands
	return CustomTask::InputParam(pName,input,gScaling);
}

// Set text-based parameters
void PeriodicXPIC::SetTextParameter(char *fxn,char *ptr)
{
	if(ptr == (char *)&gridBCOption)
	{	// bmp file name
		if(CIstrcmp(fxn,"combined")==0)
			gridBCOption = GRIDBC_COMBINED;
		else if(CIstrcmp(fxn,"velocity")==0)
			gridBCOption = GRIDBC_VELOCITY_ONLY;
		else if(CIstrcmp(fxn,"lumped")==0)
			gridBCOption = GRIDBC_LUMPED_ONLY;
		else
		{	ThrowSAXException("GridBCOption must be 'combined', 'velocity', or 'lumped'");
			return;
		}
	}
	else
		CustomTask::SetTextParameter(fxn,ptr);
}

// Called while reading when done setting all parameters
// Make and needed settings, throw SAXExecption on error
void PeriodicXPIC::Finalize(void)
{
	if(periodicXPICorder<1 && periodicFMPMorder<1)
		ThrowSAXException("PeriodicXPIC custom task must set XPIC or FMPM order to 1 or higher");
	if(periodicXPICorder>=1 && periodicFMPMorder>=1)
		ThrowSAXException("PeriodicXPIC custom task must set XPIC or FMPM order, but not both");
	
	usingFMPM = false;
	if(periodicFMPMorder>0)
	{	periodicXPICorder = periodicFMPMorder;
		usingFMPM = true;
	}
	
	// steps override times and chack for active
	bool xpicActive = false;
	
	// mechanics
	if(periodicSteps>0)
	{	periodicTime = periodicCFL = -1.;
		xpicActive = true;
	}
	else if(periodicTime>0. || periodicCFL>0.)
		xpicActive = true;

	// conduction
	if(periodicStepsConduction>0)
	{	periodicTimeConduction = periodicCFLConduction -1.;
		xpicActive = true;
	}
	else if(periodicTimeConduction>0. || periodicCFLConduction>0.)
		xpicActive = true;
	
	// diffusion
	if(periodicStepsDiffusion>0)
	{	periodicTimeDiffusion = periodicCFLDiffusion = -1.;
		xpicActive = true;
	}
	else if(periodicTimeDiffusion>0. || periodicCFLDiffusion>0.)
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

	bodyFrc.SetUsingFMPM(usingFMPM);
	if(gridBCOption<0)
		gridBCOption = usingFMPM ? GRIDBC_COMBINED : GRIDBC_LUMPED_ONLY;
	bodyFrc.SetGridBCOption(gridBCOption);
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
		if(!usingFMPM && gridBCOption==GRIDBC_VELOCITY_ONLY)
			gridBCOption = GRIDBC_LUMPED_ONLY;
		if(periodicXPICorder>0)
		{	cout << "      Grid Velocity BC Option: ";
			if(gridBCOption==GRIDBC_COMBINED)
				cout << "combined" << endl;
			else if(gridBCOption==GRIDBC_VELOCITY_ONLY)
				cout << "velocity only" << endl;
			else
				cout << "lumped only" << endl;
		}
#if MM_XPIC == 1
		if(bodyFrc.XPICVectors()>3)
			cout << "      " << GetType() << " accounts for contact" << endl;
#endif
	}

	// Conduction XPIC
	if(ConductionTask::active)
	{	cout << "   Periodic " << GetType() << " interval for conduction: ";
		if(periodicTimeConduction>0. || periodicCFLConduction>0.)
		{	if(periodicCFLConduction>0.) periodicTimeConduction = periodicCFLConduction*conduction->GetTimeStep();
			cout << periodicTimeConduction*UnitsController::Scaling(1.e3) << " " << UnitsController::Label(ALTTIME_UNITS) << endl;
			nextPeriodicTimeConduction = periodicTimeConduction;
			TransportTask::hasXPICOption = true;
			if(periodicTimeConduction<=timestep)
			{	periodicTimeConduction = -1.;
				periodicStepsConduction = 1;
				nextPeriodicStepConduction = periodicStepsConduction;
			}
		}
		else if(periodicStepsConduction>0)
		{	if(periodicStepsConduction>1)
				cout << periodicStepsConduction << " time steps" << endl;
			else
				cout << "every time step" << endl;
			nextPeriodicStepConduction = periodicStepsConduction;
			TransportTask::hasXPICOption = true;
		}
		else
			cout << "not used" << endl;
	}
	
	// Diffusion XPIC
	if(fmobj->HasFluidTransport())
	{	cout << "   Periodic " << GetType() << " interval for diffusion: ";
		if(periodicTimeDiffusion>0. || periodicCFLDiffusion>0.)
		{	if(periodicCFLDiffusion>0.) periodicTimeDiffusion = periodicCFLDiffusion*diffusion->GetTimeStep();
			cout << periodicTimeDiffusion*UnitsController::Scaling(1.e3) << " " << UnitsController::Label(ALTTIME_UNITS) << endl;
			nextPeriodicTimeDiffusion = periodicTimeDiffusion;
			TransportTask::hasXPICOption = true;
			if(periodicTimeDiffusion<=timestep)
			{	periodicTimeDiffusion = -1.;
				periodicStepsDiffusion = 1;
				nextPeriodicStepDiffusion = periodicStepsDiffusion;
			}
		}
		else if(periodicStepsDiffusion>0)
		{	if(periodicStepsDiffusion>1)
				cout << periodicStepsDiffusion << " time steps" << endl;
			else
				cout << "every time step" << endl;
			nextPeriodicStepDiffusion = periodicStepsDiffusion;
			TransportTask::hasXPICOption = true;
		}
		else
			cout << "not used" << endl;
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
	doXPICConduction = false;
	doXPICDiffusion = false;

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

	if(ConductionTask::active)
	{	if(periodicTimeConduction>0.)
		{	// controlled by step time
			if(mtime+timestep>=nextPeriodicTimeConduction)
			{	doXPICConduction = true;
				nextPeriodicTimeConduction += periodicTimeConduction;
			}
		}
		else if(periodicStepsConduction>0)
		{	// controlled by step count
			if(fmobj->mstep+1>=nextPeriodicStepConduction)
			{	doXPICConduction = true;
				nextPeriodicStepConduction += periodicStepsConduction;
			}
		}
	}
	
	if(fmobj->HasFluidTransport())
	{	if(periodicTimeDiffusion>0.)
		{	// controlled by step time
			if(mtime+timestep>=nextPeriodicTimeDiffusion)
			{	doXPICDiffusion = true;
				nextPeriodicTimeDiffusion += periodicTimeDiffusion;
			}
		}
		else if(periodicStepsDiffusion>0)
		{	// controlled by step count
			if(fmobj->mstep+1>=nextPeriodicStepDiffusion)
			{	doXPICDiffusion = true;
				nextPeriodicStepDiffusion += periodicStepsDiffusion;
			}
		}
	}
	
	return nextTask;
}

// Adjust time step now
CustomTask *PeriodicXPIC::StepCalculation(void)
{
	// reset to FLIP
	bodyFrc.SetXPICOrder(0);
	bodyFrc.SetUsingVstar(VSTAR_NOT_USED);

	TransportTask::XPICOrder = 0;
	if(ConductionTask::active) conduction->SetUsingTransportXPIC(false);
	if(fmobj->HasFluidTransport()) diffusion->SetUsingTransportXPIC(false);

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
	
	// conduction XPIC
	if(doXPICConduction)
	{	// switch to XPIC for next time step
		TransportTask::XPICOrder = periodicXPICorder;
		conduction->SetUsingTransportXPIC(true);
		
		if(verbose)
		{	cout <<"# Use " << GetType() << "(" << periodicXPICorder << ") for conduction in step " << (fmobj->mstep+1) << endl;
		}
	}
	
	// conduction XPIC
	if(doXPICDiffusion)
	{	// switch to XPIC for next time step
		TransportTask::XPICOrder = periodicXPICorder;
		diffusion->SetUsingTransportXPIC(true);
		
		if(verbose)
		{	cout <<"# Use " << GetType() << "(" << periodicXPICorder << ") for diffusion in step " << (fmobj->mstep+1) << endl;
		}
	}

	return nextTask;
}

// return the type
const char *PeriodicXPIC::GetType(void)
{   return usingFMPM ? "FMPM" : "XPIC" ;
}
