/********************************************************************************
	RunCustomTasksTask.hpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
 ********************************************************************************/

#ifndef _RUNCUSTOMTASKSTASK_

#define _RUNCUSTOMTASKSTASK_

#include "NairnMPM_Class/MPMTask.hpp"

class RunCustomTasksTask : public MPMTask
{
	public:
	
		// constructor
		RunCustomTasksTask(const char *);
	
		// required methods
		virtual void Execute(void);
	
	protected:
		double zDeriv[MaxShapeNds];
	
};

#endif

