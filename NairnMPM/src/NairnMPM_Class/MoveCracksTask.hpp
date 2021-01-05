/********************************************************************************
	MoveCracksTask.hpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _MOVECRACKSTASK_

#define _MOVECRACKSTASK_

#include "NairnMPM_Class/MPMTask.hpp"

class MoveCracksTask : public MPMTask
{
	public:
	
		// constructor
		MoveCracksTask(const char *);
	
		// required methods
		virtual bool Execute(int);
	
	protected:
	
};

#endif
