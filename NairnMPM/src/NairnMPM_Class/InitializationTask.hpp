/********************************************************************************
	InitializationTask.hpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
	MPMTask, CommonTask
********************************************************************************/

#ifndef _INITIALIZATIONTASK_

#define _INITIALIZATIONTASK_

#include "NairnMPM_Class/MPMTask.hpp"

class InitializationTask : public MPMTask
{
public:
	
	// constructor
	InitializationTask(const char *);
	
	// required methods
	virtual void Execute(void);
	
protected:
};

#endif
