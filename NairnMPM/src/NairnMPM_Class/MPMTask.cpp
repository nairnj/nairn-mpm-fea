/********************************************************************************
	MPMTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
********************************************************************************/

#include "NairnMPM_Class/MPMTask.hpp"
#include "System/ArchiveData.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Patches/GridPatch.hpp"
#include "Nodes/NodalPoint.hpp"

#pragma mark MPMTask::Constructors

// constructor
MPMTask::MPMTask(const char *name) : CommonTask(name) {}

#pragma mark MPMTask::Progress and Profiling Methods

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

// track times
void MPMTask::TrackTimes(double beginTime,double beginETime)
{	totalTaskTime += fmobj->CPUTime()-beginTime;
	totalTaskETime += fmobj->ElapsedTime()-beginETime;
}

// report on times
void MPMTask::WriteProfileResults(int nsteps,double timePerStep,double eTimePerStep)
{
	cout << "Task #" << taskNumber << ": "<< GetTaskName();
#ifdef _OPENMP
    // CPU time per step (not allows good number)
	double taskPerStep = 1000.*totalTaskTime/(double)nsteps;
	cout << ": " << taskPerStep << " CPU_ms/step (" << 100.*taskPerStep/timePerStep << "%), ";
#endif
	
    // elapsed time per step
	double eTaskPerStep = 1000.*totalTaskETime/(double)nsteps;
	cout << eTaskPerStep << " ms/step (" << 100.*eTaskPerStep/eTimePerStep << "%, "
			<< totalTaskTime/totalTaskETime << ")";
	
	cout << endl;
}

#endif

#pragma mark MPMTASK::Static Parallel Methods

// get patch number of the current thread (or 0 if not parallel
// same as thread num in loops that are not using patches
int MPMTask::GetPatchNumber(void)
{
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

// get pointer to node, which might be a ghost node in parallel code
NodalPoint *MPMTask::GetNodePointer(int pn,int nodeNum)
{
#ifdef _OPENMP
    return patches[pn]->GetNodePointer(nodeNum);
#else
    return nd[nodeNum];
#endif
}
