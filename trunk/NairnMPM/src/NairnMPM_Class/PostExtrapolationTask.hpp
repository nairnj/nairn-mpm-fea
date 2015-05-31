/********************************************************************************
	PostExtrapolationTask.hpp
	nairn-mpm-fea

	Created by John Nairn on March 6, 2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _POSTEXTRAPOLATIONTASK_

#define _POSTEXTRAPOLATIONTASK_

#include "NairnMPM_Class/MPMTask.hpp"

class PostExtrapolationTask : public MPMTask
{
	public:
	
		// constructor
		PostExtrapolationTask(const char *);
	
		// required methods
		virtual void Execute(void);
	
	protected:
};

#endif
