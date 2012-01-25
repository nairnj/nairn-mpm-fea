/********************************************************************************
	UpdateMomentaTask.hpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _UPDATEMOMENTATASK_

#define _UPDATEMOMENTATASK_

#include "NairnMPM_Class/MPMTask.hpp"

class UpdateMomentaTask : public MPMTask
{
	public:
	
		// constructor
		UpdateMomentaTask(const char *);
	
		// required methods
		virtual void Execute(void);
	
	protected:
};

#endif
