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
PeriodicXPIC::PeriodicXPIC() : CustomTask()
{
	verbose = 0;
	periodicXPICorder = 0;
	
	periodicTime = -1.;
	periodicCFL = -1.;
	periodicSteps = -1;

	periodicFMPMorder = 0;
	gridBCOption = -1;
	
	periodicTimeConduction = -1.;
	periodicCFLConduction = -1.;
    periodicStepsConduction = -1;
	conductionFractionFMPM = -1.;
	for(int i=0;i<MAX_DIFFUSION_TASKS;i++)
	{	periodicTimeDiffusion[i] = -1.;
		periodicCFLDiffusion[i] = -1.;
		diffusionFractionFMPM[i] = -1.;
        periodicStepsDiffusion[i] = -1;
	}
}

// Return name of this task
const char *PeriodicXPIC::TaskName(void) { return "Periodic XPIC Implementation"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
char *PeriodicXPIC::InputParam(char *pName,int &input,double &gScaling)
{
	// check for terminal number
	int paramNum = 0;
	char lastChar = pName[strlen(pName)-1];
	if(lastChar>='0' && lastChar<='9')
	{	paramNum = (int)lastChar - (int)'0';
		pName[strlen(pName)-1] = 0;
	}
	
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
	{	input=DOUBLE_NUM;
		return (char *)&conductionFractionFMPM;
	}
	
	// diffusion
	else if(strcmp(pName,"periodicTimeDiffusion")==0)
	{	if(paramNum>=MAX_DIFFUSION_TASKS)
		{	cout << "Diffusion parameter number must be less the maximum number of diffusion tasks";
			cout << " (" << MAX_DIFFUSION_TASKS << ")" << endl;
			return NULL;
		}
		input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&periodicTimeDiffusion[paramNum],gScaling,1.e-3);
	}
	
	else if(strcmp(pName,"periodicCFLDiffusion")==0)
	{	if(paramNum>=MAX_DIFFUSION_TASKS)
		{	cout << "Diffusion parameter number must be less the maximum number of diffusion tasks";
			cout << " (" << MAX_DIFFUSION_TASKS << ")" << endl;
			return NULL;
		}
		input=DOUBLE_NUM;
		return (char *)&periodicCFLDiffusion[paramNum];
	}
	
	else if(strcmp(pName,"periodicStepsDiffusion")==0)
	{	if(paramNum>=MAX_DIFFUSION_TASKS)
		{	cout << "Diffusion parameter number must be less the maximum number of diffusion tasks";
			cout << " (" << MAX_DIFFUSION_TASKS << ")" << endl;
			return NULL;
		}
		input=DOUBLE_NUM;
		return (char *)&diffusionFractionFMPM[paramNum];
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
			ThrowSAXException("GridBCOption must be 'combined', 'velocity', or 'lumped'");
	}
	else
		CustomTask::SetTextParameter(fxn,ptr);
}

// Called while reading when done setting all parameters
// Make any needed settings, throw SAXExecption on error
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
	if(conductionFractionFMPM>0.)
	{	periodicTimeConduction = periodicCFLConduction = -1.;
		if(conductionFractionFMPM<1.)
		{	periodicStepsConduction=1;
		}
		else
		{	periodicStepsConduction=int(conductionFractionFMPM+0.5);
			conductionFractionFMPM=1.;
		}
		xpicActive = true;
	}
	else if(periodicTimeConduction>0. || periodicCFLConduction>0.)
		xpicActive = true;
	
	// diffusion tasks
	for(int i=0;i<MAX_DIFFUSION_TASKS;i++)
    {   if(diffusionFractionFMPM[i]>0.)
		{	periodicTimeDiffusion[i] = periodicCFLDiffusion[i] = -1.;
			if(diffusionFractionFMPM[i]<1.)
			{	periodicStepsDiffusion[i]=1;
			}
			else
			{	periodicStepsDiffusion[i]=int(diffusionFractionFMPM[i]+0.5);
				diffusionFractionFMPM[i]=1.;
			}
			xpicActive = true;
		}
		else if(periodicTimeDiffusion[i]>0. || periodicCFLDiffusion[i]>0.)
        {   periodicStepsDiffusion[i] = -1;
            diffusionFractionFMPM[i] = 1.;          // need because not used in this mode
			xpicActive = true;
        }
	}
	
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
        // if this time < timestep, then done every time step
		nextPeriodicTime = periodicTime;
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
            // if this time < timestep, then done every time step
            nextPeriodicTimeConduction = periodicTimeConduction;
			TransportTask::hasXPICOption = true;
		}
		else if(periodicStepsConduction>0)
		{	if(periodicStepsConduction>1)
				cout << periodicStepsConduction << " time steps" << endl;
			else if(conductionFractionFMPM<1.)
				cout << "fraction " << conductionFractionFMPM << " with rest FLIP every time step" << endl;
			else
				cout << "every time step" << endl;
			nextPeriodicStepConduction = periodicStepsConduction;
			TransportTask::hasXPICOption = true;
		}
		else
			cout << "not used" << endl;
	}
	
	// Diffusion XPIC/FMPM (each can set their own)
	// Note: assumes diffusion tasks created before creating this custom task
	// [0] is always for standard diffusion (if active)
	// [i>0] is for otherDiffsion tasks (numbered from 1)
	if(numDiffusion>0)
    {	DiffusionTask *nextTask;
        int d1,d2;
        if(diffusion==NULL)
        {   nextTask = otherDiffusion;
            d1 = 1;
            d2 = numDiffusion+1;
        }
        else
        {   nextTask = diffusion;
            d1 = 0;
            d2 = numDiffusion;
        }
		for(int i=d1;i<d2;i++)
		{	cout << "   Periodic " << GetType() << " interval for ";
			cout << nextTask->StyleName() << " (number " << nextTask->GetNumber() << "): ";
			if(periodicTimeDiffusion[i]>0. || periodicCFLDiffusion[i]>0.)
			{	if(periodicCFLDiffusion[i]>0.)
					periodicTimeDiffusion[i] = periodicCFLDiffusion[i]*nextTask->GetTimeStep();
				cout << periodicTimeDiffusion[i]*UnitsController::Scaling(1.e3) << " "
							<< UnitsController::Label(ALTTIME_UNITS) << endl;
                // if this time < timestep, then done every time step
				nextPeriodicTimeDiffusion[i] = periodicTimeDiffusion[i];
                TransportTask::hasXPICOption = true;
			}
			else if(periodicStepsDiffusion[i]>0)
			{	if(periodicStepsDiffusion[i]>1)
					cout << periodicStepsDiffusion[i] << " time steps" << endl;
				else if(diffusionFractionFMPM[i]<1.)
					cout << "fraction " << diffusionFractionFMPM[i] << " with rest FLIP every time step" << endl;
				else
					cout << "every time step" << endl;
				nextPeriodicStepDiffusion[i] = periodicStepsDiffusion[i];
				TransportTask::hasXPICOption = true;
			}
			else
				cout << "not used" << endl;
			
			// next task
			nextTask = i==0 ? otherDiffusion : (DiffusionTask *)nextTask->GetNextTransportTask();
		}
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
	for(int i=0;i<numDiffusion;i++) doXPICDiffusion[i] = false;

	if(periodicTime>0.)
	{	// controlledby step time
		if(mtime+timestep>=nextPeriodicTime)
		{	doXPIC = true;
			nextPeriodicTime = mtime+timestep+periodicTime;
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
				nextPeriodicTimeConduction = mtime+timestep+periodicTimeConduction;
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
	
	// set true if any diffusion task wants XPIC.FMPM
	if(numDiffusion>0)
	{	DiffusionTask *nextTask;
        int d1,d2;
        if(diffusion==NULL)
        {   nextTask = otherDiffusion;
            d1 = 1;
            d2 = numDiffusion+1;
        }
        else
        {   nextTask = diffusion;
            d1 = 0;
            d2 = numDiffusion;
        }
        for(int i=d1;i<d2;i++)
		{	if(periodicTimeDiffusion[i]>0.)
			{	// controlled by step time
				if(mtime+timestep>=nextPeriodicTimeDiffusion[i])
				{	doXPICDiffusion[i] = true;
					nextPeriodicTimeDiffusion[i] = mtime+timestep+periodicTimeDiffusion[i];
				}
			}
			else if(periodicStepsDiffusion[i]>0)
			{	// controlled by step count
				if(fmobj->mstep+1>=nextPeriodicStepDiffusion[i])
				{	doXPICDiffusion[i] = true;
					nextPeriodicStepDiffusion[i] += periodicStepsDiffusion[i];
				}
			}
			
			nextTask = (i==0) ? otherDiffusion : (DiffusionTask *)nextTask->GetNextTransportTask();
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
	if(ConductionTask::active) conduction->SetUsingTransportXPIC(false,1.);
	if(numDiffusion>0) DiffusionTask::SetDiffusionXPIC(false);

	// mechanics XPIC/FMPM
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
		{	// for order 1
			bodyFrc.SetUsingVstar(VSTAR_NOT_USED);
		}
	
		if(verbose)
		{	cout <<"# Use " << GetType() << "(" << periodicXPICorder << ") in step " << (fmobj->mstep+1) << endl;
		}
	}
	
	// conduction XPIC
	if(doXPICConduction)
	{	// switch to XPIC for next time step
		TransportTask::XPICOrder = periodicXPICorder;
		conduction->SetUsingTransportXPIC(true,conductionFractionFMPM);
		
		if(verbose)
		{	cout <<"# Use " << GetType() << "(" << periodicXPICorder << ") for conduction in step " << (fmobj->mstep+1) << endl;
		}
	}
	
	// diffusion XPIC/FMPM
	if(numDiffusion>0)
	{	DiffusionTask *nextTask;
        int d1,d2;
        if(diffusion==NULL)
        {   nextTask = otherDiffusion;
            d1 = 1;
            d2 = numDiffusion+1;
        }
        else
        {   nextTask = diffusion;
            d1 = 0;
            d2 = numDiffusion;
        }
        for(int i=d1;i<d2;i++)
		{	if(doXPICDiffusion[i])
			{	// switch to XPIC for next time step (make sure fraction is 1 unless mean it)
				TransportTask::XPICOrder = periodicXPICorder;
				nextTask->SetUsingTransportXPIC(true,diffusionFractionFMPM[i]);
				
				if(verbose)
				{	cout <<"# Use " << GetType() << "(" << periodicXPICorder << ") for ";
					cout << nextTask->TaskName() << " in step " << (fmobj->mstep+1) << endl;
				}
			}
			
			nextTask = (i==0) ? otherDiffusion : (DiffusionTask *)nextTask->GetNextTransportTask();
		}
	}

	return nextTask;
}

// return the type
const char *PeriodicXPIC::GetType(void)
{   return usingFMPM ? "FMPM" : "XPIC" ;
}
