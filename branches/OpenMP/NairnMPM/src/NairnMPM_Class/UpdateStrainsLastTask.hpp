/********************************************************************************
	UpdateStrainsLastTask.hpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _UPDATESTRAINSLASTTASK_

#define _UPDATESTRAINSLASTTASK_

#include "NairnMPM_Class/MPMTask.hpp"

class UpdateStrainsLastTask : public MPMTask
{
	public:
	
		// constructor
		UpdateStrainsLastTask(const char *);
	
		// required methods
		virtual void Execute(void);
	
	protected:
		double zDeriv[MaxShapeNds];
	
};

#endif
