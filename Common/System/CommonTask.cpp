/********************************************************************************
	CommonTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.

********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "CommonTask.hpp"
#include "System/MPMPrefix.hpp"

// class statics
int CommonTask::numberOfTasks=0;

#pragma mark CONSTRUCTORS

// empty constructor
CommonTask::CommonTask() {}

// default constructor
CommonTask::CommonTask(const char *name)
{
	numberOfTasks++;
	taskNumber=numberOfTasks;
	
	taskName=NULL;
	SetTaskName(name);
	
	nextTask=NULL;
	
	totalTaskTime=0.;
	totalTaskETime=0.;
}

CommonTask::~CommonTask()
{
	if(taskName!=NULL) delete [] taskName;
}

#pragma mark ACCESSORS

// Set task name
// throws std::bad_alloc
void CommonTask::SetTaskName(const char *descrip)
{	if(taskName!=NULL) delete [] taskName;
	if(descrip==NULL)
	{	taskName=NULL;
		return;
	}
    taskName=new char[strlen(descrip)+1];
    strcpy(taskName,descrip);
}
char *CommonTask::GetTaskName(void) { return taskName; }

// next task
CommonTask *CommonTask::GetNextTask(void) { return nextTask; }
void CommonTask::SetNextTask(CommonTask *obj) { nextTask=obj; }

