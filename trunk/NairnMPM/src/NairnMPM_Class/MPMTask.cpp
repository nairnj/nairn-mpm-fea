/********************************************************************************
	MPMTask.cpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
********************************************************************************/

#include "NairnMPM_Class/MPMTask.hpp"
#include "System/ArchiveData.hpp"

#pragma mark CONSTRUCTORS

// constructor
MPMTask::MPMTask(const char *name) : CommonTask(name) {}

#ifdef LOG_PROGRESS

// for debugging
void MPMTask::WriteLogFile(void)
{	// task number
	char logLine[200];
	sprintf(logLine,"  Task #%d: %s",taskNumber,GetTaskName());
	archiver->WriteLogFile(logLine,NULL);
}

#endif

#ifdef _PROFILE_TASKS_

// report on times
void MPMTask::WriteProfileResults(int nsteps,double timePerStep)
{
	cout << "Task #" << taskNumber << ": "<< GetTaskName();
	double taskPerStep=1000.*totalTaskTime/(double)nsteps;
	cout << ": " << taskPerStep << " ms/step (" << 100.*taskPerStep/timePerStep << "%)" << endl;
}

#endif

