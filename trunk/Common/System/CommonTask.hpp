/********************************************************************************
	CommonTask.hpp
	NairnFEA and NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
		none
********************************************************************************/

#ifndef _COMMONTASK_

#define _COMMONTASK_

class CommonTask
{
	public:
	
		//  Constructors and Destructor
		CommonTask();
		CommonTask(const char *);
		virtual ~CommonTask();
	
		// abstract methods
		virtual void Execute(void) = 0;
	
		// accessors
		void SetTaskName(const char *);
		char *GetTaskName(void);
		void SetNextTask(CommonTask *);
		CommonTask *GetNextTask(void);
	
	protected:
		char *taskName;
		CommonTask *nextTask;
		static int numberOfTasks;
		int taskNumber;
		double totalTaskTime;
};

#endif