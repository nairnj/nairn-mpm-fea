/******************************************************************************** 
    CustomTask.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Fri Aug 15 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Custom_Tasks/CalcJKTask.hpp"

// globals
CustomTask *theTasks=NULL;
int CustomTask::numberCustomHistoryVariables = 0;

#pragma mark INITIALIZE

// Constructors
CustomTask::CustomTask() { nextTask=NULL; }

// Destructor (and it is virtual)
CustomTask::~CustomTask() { }

// Return name of this task
const char *CustomTask::TaskName(void) { return "Unnamed Custom Task"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
// On error, cout a message and return NULL or throw SAXException()
char *CustomTask::InputParam(char *pName,int &input,double &gScaling) { return NULL; }

// If pointer type in InputParam is set to TEXT_PARAMETER, this method is called
//      immediately after with the text of the parameter
// If value not accepted, throw exception
// throws SAXException()
void CustomTask::SetTextParameter(char *value,char *ptr)
{	
	ThrowSAXException("The text parameter with value '%s' is not valid for current custom task",value);
}

#pragma mark GENERIC TASK METHODS

// Called while reading when done setting all parameters
// Make and needed settings, throw SAXExecption on error
void CustomTask::Finalize(void) {}

// called once at start of MPM analysis in CreateTasks (which is after Preliminary calcs)
// Print details on the tasks (errors might be better flagged in Finalize()
CustomTask *CustomTask::Initialize(void)
{
    cout << "Unknown custom task" << endl;
    return nextTask;
}

// Increment history by number this task needs
// Also store which ones this task will be using
void CustomTask::ClaimHistoryVariables() {}

// called when MPM step is getting ready to do custom tasks
// return nextTask when done
CustomTask *CustomTask::PrepareForStep(bool &doNodalExtraps)
{ return nextTask; }

// do custom calculation
CustomTask *CustomTask::StepCalculation(void) { return nextTask; }

// do custom calculations before the first time step
CustomTask *CustomTask::Step0Calculation(void) { return nextTask; }

// Called when custom tasks are all done on a step
CustomTask *CustomTask::FinishForStep(bool &removeMe) { return nextTask; }

// called after the analysis is done. A task can output a report to the main.mpm
// file. It should end with an empty line. Task with a report must override
// both methods and return true is HasReport()
bool CustomTask::HasReport(void) { return false; }
CustomTask *CustomTask::Report(void) { return nextTask; }

#pragma mark TASK EXTRAPOLATION METHODS

// initialize prior to node and particle extrapolations
CustomTask *CustomTask::BeginExtrapolations(void) { return nextTask; }

// called when done with node and particle extrapolations
CustomTask *CustomTask::EndExtrapolations(void) { return nextTask; }

// add particle data to a node during extrapolations (wt is mp*fn[i])
CustomTask *CustomTask::NodalExtrapolation(NodalPoint *ndmi,MPMBase *mpnt,short vfld,int matfld,double wt,short isRigid)
{ return nextTask; }



