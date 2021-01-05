/********************************************************************************
	ExtrapolateRigidBCsTask.hpp
	nairn-mpm-fea

	Created by John Nairn on March 6, 2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _EXTRAPOLATERIGIDBCSTASK_

#define _EXTRAPOLATERIGIDBCSTASK_

class BoundaryCondition;

#include "NairnMPM_Class/MPMTask.hpp"

class ExtrapolateRigidBCsTask : public MPMTask
{
	public:
	
		// constructor
		ExtrapolateRigidBCsTask(const char *);
	
		// required methods
		virtual bool Execute(int);
	
	protected:
};

#endif
