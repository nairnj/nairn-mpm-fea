/********************************************************************************
	UpdateStrainsFirstTask.hpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _UPDATESTRAINSFIRSTTASK_

#define _UPDATESTRAINSFIRSTTASK_

#include "NairnMPM_Class/MPMTask.hpp"

class UpdateStrainsFirstTask : public MPMTask
{
	public:
	
		// constructor
		UpdateStrainsFirstTask(const char *);
	
		// required methods
		virtual void Execute(void);
	
        // class methods
        static void FullStrainUpdate(double,int,int);
	
	protected:
};

#endif
