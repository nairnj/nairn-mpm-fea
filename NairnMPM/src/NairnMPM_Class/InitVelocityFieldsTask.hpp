/********************************************************************************
	InitVelocityFieldsTask.hpp
	nairn-mpm-fea

	Created by John Nairn on March 5, 2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _INITVELOCITYFIELDSTASK_

#define _INITVELOCITYFIELDSTASK_

#include "NairnMPM_Class/MPMTask.hpp"

class InitVelocityFieldsTask : public MPMTask
{
	public:
#ifdef PREHASH_CRACKS
        static bool prehashed;
#endif
	
		// constructor
		InitVelocityFieldsTask(const char *);
	
		// required methods
		virtual bool Execute(int);
	
	protected:
};

#endif
