/********************************************************************************
    AdjustTimeStepTask.cpp
    nairn-mpm-fea

    Created by John Nairn on 9/24/12.
    Copyright (c) 2012 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Custom_Tasks/AdjustTimeStepTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "System/ArchiveData.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "Elements/ElementBase.hpp"
#include "System/UnitsController.hpp"

#pragma mark Constructors and Destructors

// Constructors
AdjustTimeStepTask::AdjustTimeStepTask()
{
	customAdjustTime = -1.;
	nextCustomAdjustTime = -1.;
    verbose = 0;
    lastReportedTimeStep = -1;
	velocityCFL = -1.;
}

// Return name of this task
const char *AdjustTimeStepTask::TaskName(void) { return "Periodically adjust MPM time step"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
char *AdjustTimeStepTask::InputParam(char *pName,int &input,double &gScaling)
{
    if(strcmp(pName,"adjustTime")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&customAdjustTime,gScaling,1.e-3);
    }
		
    else if(strcmp(pName,"verbose")==0)
    {	input=INT_NUM;
        return (char *)&verbose;
    }
	
    else if(strcmp(pName,"velocityCFL")==0)
    {	input=DOUBLE_NUM;
        return (char *)&velocityCFL;
    }
	
	// check remaining commands
    return CustomTask::InputParam(pName,input,gScaling);
}

#pragma mark GENERIC TASK METHODS

// called once at start of MPM analysis - initialize and print info
CustomTask *AdjustTimeStepTask::Initialize(void)
{
    cout << "Periodically adjust MPM time step." << endl;
	
	// time interval
	cout << "   Adjust interval: ";
	if(customAdjustTime>=0.)
	{	cout << customAdjustTime*UnitsController::Scaling(1.e3) << " " << UnitsController::Label(ALTTIME_UNITS) << endl;
		nextCustomAdjustTime = customAdjustTime;
	}
	else
		cout << "same as particle archives" << endl;
	
	// CFL factor - output and change to relative value
	if(velocityCFL>0.)
	{	cout << "   CFL factor on particle velocities: " << velocityCFL << endl;
		velocityCFL /= fmobj->GetCFLCondition();
	}
	else
		velocityCFL = 1.;
	
    if(verbose!=0)
        cout << "   Verbose: yes" << endl;
    else
        cout << "   Verbose: no" << endl;
	
    return nextTask;
}

// called when MPM step is getting ready to do custom tasks
// never uses extrapolations so no need to set needExtraps
CustomTask *AdjustTimeStepTask::PrepareForStep(bool &needExtraps)
{
	if(customAdjustTime>=0.)
	{	if(mtime+timestep>=nextCustomAdjustTime)
        {	doAdjust = true;
            nextCustomAdjustTime += customAdjustTime;
        }
        else
            doAdjust = false;
	}
	else
		doAdjust = archiver->WillArchive();
    
    return nextTask;
}

// Adjust time step now
CustomTask *AdjustTimeStepTask::StepCalculation(void)
{
    // exit when not needed
    if(!doAdjust) return nextTask;
    
    // get grid dimensions
    double dcell = mpmgrid.GetMinCellDimension();
    
    if(lastReportedTimeStep<0) lastReportedTimeStep = timestep;
    
    // set globals
    bool Max_Velocity_Condition = false;
    timestep = 1.e15;
    propTime = 1.e15;
	
     // loop over nonrigid material points
    for(int p=0;p<nmpmsNR;p++)
	{	// material id
		short matid=mpm[p]->MatID();
        
        // check time step using convergence condition (wave speed of material)
        double crot = theMaterials[matid]->CurrentWaveSpeed(fmobj->IsThreeD(),mpm[p],0);
		
		// check to see if particle speed is above critical speed
		double crot1 = fabs(mpm[p]->vel.x)+fabs(mpm[p]->vel.y);
		if(fmobj->IsThreeD()) crot1 += fabs(mpm[p]->vel.z);
		crot1 /= velocityCFL;
		
		// Pick highest speed
		if(crot1>crot) {
			crot=crot1;
			Max_Velocity_Condition = true;
		}
		
		// test time step
		double tst = fmobj->GetCFLCondition()*dcell/crot;
        if(tst<timestep) timestep = tst;
        
        // propagation time (in sec)
        tst = fmobj->GetPropagationCFLCondition()*dcell/crot;
        if(tst<propTime) propTime = tst;
	}
	
    // verify time step and make smaller if needed
	strainTimestep = (fmobj->mpmApproach==USAVG_METHOD) ? timestep/2. : timestep ;
    
	// propagation time step (no less than timestep)
    if(propTime<timestep) propTime = timestep;
    
    // report if changed by 5% since last reported change
    if(verbose!=0)
    {   double ratio = timestep/lastReportedTimeStep;
        if(ratio < 0.95 || ratio>1.05)
		{	if(timestep<lastReportedTimeStep)
				cout << "# time step reduced to " << timestep*UnitsController::Scaling(1.e3) << " " << UnitsController::Label(ALTTIME_UNITS);
			else
				cout << "# time step increased to " << timestep*UnitsController::Scaling(1.e3) << " " << UnitsController::Label(ALTTIME_UNITS);
			if(Max_Velocity_Condition)
				cout << " (some velocities exceed wave speed)";
			cout << endl;
		}
        lastReportedTimeStep = timestep;
    }

    return nextTask;
}

