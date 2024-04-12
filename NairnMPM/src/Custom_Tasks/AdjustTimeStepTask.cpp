/********************************************************************************
    AdjustTimeStepTask.cpp
    nairn-mpm-fea

    Created by John Nairn on 9/24/12.
    Copyright (c) 2012 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Custom_Tasks/AdjustTimeStepTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "System/ArchiveData.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "Elements/ElementBase.hpp"
#include "System/UnitsController.hpp"
#include "Exceptions/CommonException.hpp"

// save base tranport time step
double AdjustTimeStepTask::transportBaseTimeStep = 1.e30;

// class globals
AdjustTimeStepTask *adjustTimeStepTask=NULL;

double minTimeStep=1.e15,maxTimeStep=-1.;

#pragma mark Constructors and Destructors

// Constructors
AdjustTimeStepTask::AdjustTimeStepTask() : CustomTask()
{
	customAdjustTime = -1.;
	nextCustomAdjustTime = -1.;
    verbose = 0;
    lastReportedTimeStep = -1;
	velocityCFL = -1.;
    reportRatio = 0.;
    maxIncrease = -1.;
    checkTransportTimeStep = 0;   // normally not needed
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
    {	input=DOUBLE_NUM;
        return (char *)&reportRatio;
    }
	
    else if(strcmp(pName,"velocityCFL")==0)
    {	input=DOUBLE_NUM;
        return (char *)&velocityCFL;
    }
	
    else if(strcmp(pName,"maxIncrease")==0)
    {   input=DOUBLE_NUM;
        return (char *)&maxIncrease;
    }
    
    else if(strcmp(pName,"checkTransportTimeStep")==0)
    {   input=INT_NUM;
        return (char *)&checkTransportTimeStep;
    }
    
	// check remaining commands
    return CustomTask::InputParam(pName,input,gScaling);
}

#pragma mark GENERIC TASK METHODS

// called once at start of MPM analysis - initialize and print info
CustomTask *AdjustTimeStepTask::Initialize(void)
{
    // save task to a global
    if(adjustTimeStepTask!=NULL)
        throw CommonException("Only one AdjustTimeStepTask allowed in a simulations","AdjustTimeStepTask::Initialize");
    adjustTimeStepTask = this;

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
    
    if(maxIncrease>0.)
        cout << "   Time step increases limited to " << (100.*maxIncrease) << "%" << endl;
	
    // set verbose from entered reportRatio
    // To support all modes, 0 is not verbose, 1 is verbose with reportRatio 2
    if(reportRatio<=0.99)
        verbose = 0;
    else
    {	verbose = 1;
        if(reportRatio<1.01) reportRatio = 2.;
    }
    
    if(verbose!=0)
        cout << "   Verbose: yes (reporting " << reportRatio << "X changes)" << endl;
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
    
    // if not reported yet, choose current time step
    if(lastReportedTimeStep<0)
        SetLastReportedTimeStep(timestep);
    
    // set globals
    double maxTimeStep=1.e15,maxPropTime=1.e15;
    if(maxIncrease>0.)
    {   maxTimeStep = (1.+maxIncrease)*timestep;
        maxPropTime = (1.+maxIncrease)*propTime;
    }
	double dcell = mpmgrid.GetMinCellDimension();
    int Max_Velocity_Condition = 0;
    double newTimestep = 1.e15;
    double newPropTime = 1.e15;
	
    // loop over nonrigid material points
#pragma omp for
   for(int p=0;p<nmpmsNR;p++)
   {   if(mpm[p]->InReservoir()) continue;
       // material id
       short matid=mpm[p]->MatID();

#ifndef TRANSPORT_ONLY
       // check time step using convergence condition (wave speed of material)
       double crot = theMaterials[matid]->CurrentWaveSpeed(fmobj->IsThreeD(),mpm[p],0);
       
       // check to see if particle speed is above critical speed
       double crot1 = fabs(mpm[p]->vel.x)+fabs(mpm[p]->vel.y);
       if(fmobj->IsThreeD()) crot1 += fabs(mpm[p]->vel.z);
       crot1 /= velocityCFL;
       
       // Pick highest speed
       if(crot1>crot)
       {   crot=crot1;
#pragma omp atomic
           Max_Velocity_Condition++;
       }
       
       // test time step
       double tst = fmin(fmobj->GetCFLCondition()*dcell/crot,maxTimeStep);
       if(tst<newTimestep)
       {
#pragma omp critical (adjustimestep)
           {
               newTimestep = tst;
           }
       }
       
       // propagation time (in sec)
       tst = fmin(fmobj->GetPropagationCFLCondition()*dcell/crot,maxPropTime);
       if(tst<newPropTime)
       {
#pragma omp critical (adjustimestep)
           {
               newPropTime = tst;
           }
       }
#endif
       
       // transport time steps (if it can change during calculations)
       if(checkTransportTimeStep!=0)
       {    // future code to check transport time steps
       }
    }
    
    // in case no particles checked for time step
    if(newTimestep>1.e14) return nextTask;
    
    // If each particle not checked above for transport, make duer
    // does not exceed the base transport time steo found in Preliminary calcs
    if(checkTransportTimeStep==0)
    {   if(transportBaseTimeStep<newTimestep)
            newTimestep = transportBaseTimeStep;
    }
    
    // change to new values
    ChangeTimestep(newTimestep,newPropTime,true);
	
    // report if changed by reportRatio since last reported change
    if(verbose!=0)
    {   double ratio = timestep/lastReportedTimeStep;
        if(ratio>reportRatio || ratio<1./reportRatio)
        {    cout << "# Step: " << fmobj->mstep << ": time step changed " << ratio << "X to "
                    << timestep*UnitsController::Scaling(1.e3) << " "
                    << UnitsController::Label(ALTTIME_UNITS);
            if(Max_Velocity_Condition>0)
                cout << " (some velocities exceed wave speed)";
            cout << endl;
            SetLastReportedTimeStep(timestep);
        }
    }

    return nextTask;
}

// set last reported time step
void AdjustTimeStepTask::SetLastReportedTimeStep(double newReportTime) { lastReportedTimeStep = newReportTime; }

// called after the analysis is done. A task can output a report to the main.mpm
// file. It should end with an empty line
bool AdjustTimeStepTask::HasReport(void) { return true; }
CustomTask *AdjustTimeStepTask::Report(void)
{
    // exit if never changed
    if(maxTimeStep<0.) return nextTask;
    
    cout << TaskName() << endl;
    cout << "    Time step was varied from: " << minTimeStep*UnitsController::Scaling(1.e3) << " to "
        << maxTimeStep*UnitsController::Scaling(1.e3) << " "
        << UnitsController::Label(ALTTIME_UNITS) << endl;
    cout << "    Max/Min range = " << (maxTimeStep/minTimeStep) << endl;
    cout << endl;
    return nextTask;
}

// class method to change time step from anywhere in code
void AdjustTimeStepTask::ChangeTimestep(double newTimestep,double newPropTime,bool inAdjustTask)
{
    // change time steps
    timestep = newTimestep;
    propTime = newPropTime;

#ifndef TRANSPORT_ONLY
    // verify time step and make smaller if needed
    if(fmobj->mpmApproach==USAVG_METHOD)
    {   strainTimestepFirst = fractionUSF*timestep;
        strainTimestepLast = timestep - strainTimestepFirst;
    }
    else
    {   strainTimestepFirst = timestep;
        strainTimestepLast = timestep;
    }
    
    // propagation time step (no less than timestep)
    if(propTime<timestep) propTime = timestep;
#endif
    
    // tell task that outside code made a change
    if(!inAdjustTask && adjustTimeStepTask!=NULL)
        adjustTimeStepTask->SetLastReportedTimeStep(timestep);
    
    // save limits
    minTimeStep = fmin(timestep,minTimeStep);
    maxTimeStep = fmax(timestep,maxTimeStep);
}
