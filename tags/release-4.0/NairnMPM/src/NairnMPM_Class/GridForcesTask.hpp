/********************************************************************************
	GridForcesTask.hpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _GRIDFORCESTASK_

#define _GRIDFORCESTASK_

#include "NairnMPM_Class/MPMTask.hpp"

class GridForcesTask : public MPMTask
{
	public:
	
		// constructor
		GridForcesTask(const char *);
	
		// required methods
		virtual void Execute(void);
	
	protected:
		double zDeriv[MaxShapeNds];
	
};

#endif
