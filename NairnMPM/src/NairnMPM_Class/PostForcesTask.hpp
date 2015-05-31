/********************************************************************************
	PostForcesTask.hpp
	nairn-mpm-fea

	Created by John Nairn on March 8, 2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _POSTFORCESTASK_

#define _POSTFORCESTASK_

#include "NairnMPM_Class/MPMTask.hpp"

class PostForcesTask : public MPMTask
{
	public:
	
		// constructor
		PostForcesTask(const char *);
	
		// required methods
		virtual void Execute(void);
	
	protected:
};

#endif
