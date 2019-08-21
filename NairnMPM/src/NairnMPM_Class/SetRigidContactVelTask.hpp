/********************************************************************************
	SetReigidContactVelTask.hpp
	nairn-mpm-fea

	Created by John Nairn on May 12, 2016
	Copyright (c) 2016 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _SETRIGIDCONTACTVELTASK_

#define _SETRIGIDCONTACTVELTASK_

#include "NairnMPM_Class/MPMTask.hpp"

class SetRigidContactVelTask : public MPMTask
{
	public:
	
		// constructor
		SetRigidContactVelTask(const char *);
	
		// required methods
		virtual void Execute(int);
	
	protected:
};

#endif
