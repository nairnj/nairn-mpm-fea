/******************************************************************************** 
    TransportTask.cpp
    NairnMPM
    
    Created by John Nairn on 7/18/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Custom_Tasks/TransportTask.hpp"

TransportTask *transportTasks=NULL;

#pragma mark INITIALIZE

// Constructors
TransportTask::TransportTask()
{	nextTask=NULL;
}

// Destructor (and it is virtual)
TransportTask::~TransportTask() { }

#pragma mark STANDARD METHODS

// transport analysis settings
TransportTask *TransportTask::TransportOutput(void)
{	cout << "Coupled " << TaskName() << endl;
	return nextTask;
}
