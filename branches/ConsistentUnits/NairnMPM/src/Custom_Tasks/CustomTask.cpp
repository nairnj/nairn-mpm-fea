/******************************************************************************** 
    CustomTask.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Fri Aug 15 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Custom_Tasks/CalcJKTask.hpp"

// globals
CustomTask *theTasks=NULL;

#pragma mark INITIALIZE

// Constructors
CustomTask::CustomTask() { nextTask=NULL; }

// Destructor (and it is virtual)
CustomTask::~CustomTask() { }

// Return name of this task
const char *CustomTask::TaskName(void) { return "Unnamed Custom Task"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
char *CustomTask::InputParam(char *pName,int &input,double &gScaling) { return NULL; }

#pragma mark GENERIC TASK METHODS

// called once at start of MPM analysis - initialize and print info
CustomTask *CustomTask::Initialize(void)
{
    cout << "Unknown custom task" << endl;
    return nextTask;
}

// called when MPM step is getting ready to do custom tasks
// return nextTask when done
CustomTask *CustomTask::PrepareForStep(bool &doNodalExtraps)
{ return nextTask; }

// do custom calculation
CustomTask *CustomTask::StepCalculation(void) { return nextTask; }

// Called when custom tasks are all done on a step
CustomTask *CustomTask::FinishForStep(void) { return nextTask; }

#pragma mark TASK EXTRAPOLATION METHODS

// initialize prior to node and particle extrapolations
CustomTask *CustomTask::BeginExtrapolations(void) { return nextTask; }

// called when done with node and particle extrapolations
CustomTask *CustomTask::EndExtrapolations(void) { return nextTask; }

// add particle data to a node during extrapolations (wt is mp*fn[i])
CustomTask *CustomTask::NodalExtrapolation(NodalPoint *ndmi,MPMBase *mpnt,short vfld,int matfld,double wt,short isRigid)
{ return nextTask; }



