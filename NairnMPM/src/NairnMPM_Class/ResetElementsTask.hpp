/********************************************************************************
	ResetElementsTask.hpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _RESETELEMENTSTASK_

#define _RESETELEMENTSTASK_

#include "NairnMPM_Class/MPMTask.hpp"

class MPMBase;

class ResetElementsTask : public MPMTask
{
	public:
	
		// constructor
		ResetElementsTask(const char *);
	
		// required methods
		virtual void Execute(void);
	
		// class method
		static int ResetElement(MPMBase *);
	
	protected:
		void ReturnToElement(MPMBase *);

};

#endif
