/********************************************************************************
	 MPMTask.hpp
	 NairnMPM
	 
	 Created by John Nairn on July 22, 2010
	 Copyright (c) 2010 John A. Nairn, All rights reserved.
	 
	 Dependencies
		CommonTask
********************************************************************************/

#ifndef _MPMTASK_

#define _MPMTASK_

// activate to track time in each task and print results at the end
//#define _PROFILE_TASKS_

#include "System/CommonTask.hpp"

class MPMTask : public CommonTask
{
	public:
	
		// constructor
		MPMTask(const char *);
	
		// methods
#ifdef LOG_PROGRESS
		void WriteLogFile(void);
#endif
	
#ifdef _PROFILE_TASKS_
		void WriteProfileResults(int,double);
#endif
	
	protected:
	
};

#endif