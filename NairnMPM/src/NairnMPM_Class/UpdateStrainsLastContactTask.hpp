/********************************************************************************
	UpdateStrainsLastContactTask.hpp
	nairn-mpm-fea

	Created by John Nairn on 4/7/15
	Copyright (c) 2015 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _UPDATESTRAINSLASTCONTACTTASK_

#define _UPDATESTRAINSLASTCONTACTTASK_

#include "NairnMPM_Class/MassAndMomentumTask.hpp"

class UpdateStrainsLastContactTask : public MassAndMomentumTask
{
	public:
	
		// constructor
		UpdateStrainsLastContactTask(const char *);
	
		// required methods
		virtual void Execute(void);
	
};

#endif
