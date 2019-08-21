/********************************************************************************
	UpdateParticlesTask.hpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _UPDATEPARTICLESTASK_

#define _UPDATEPARTICLESTASK_

#include "NairnMPM_Class/MPMTask.hpp"

class UpdateParticlesTask : public MPMTask
{
	public:
	
		// constructor
		UpdateParticlesTask(const char *);
	
		// required methods
		virtual void Execute(int);
	
	protected:
	
};

#endif
