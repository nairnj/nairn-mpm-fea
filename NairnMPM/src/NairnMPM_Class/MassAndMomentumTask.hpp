/********************************************************************************
	MassAndMomentumTask.hpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _MASSANDMOMENTUMTASK_

#define _MASSANDMOMENTUMTASK_

class BoundaryCondition;
class MPMBase;

#include "NairnMPM_Class/MPMTask.hpp"

class MassAndMomentumTask : public MPMTask
{
	public:
	
		// constructor
		MassAndMomentumTask(const char *);
	
		// required methods
		virtual void Execute(void);
		int GetParticleFunctions(MPMBase *,int **,double *,double *,double *,double *);
	
	protected:
};

#endif
