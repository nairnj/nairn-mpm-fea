/********************************************************************************
	MassAndMomentumTask.hpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _MASSANDMOMENTUMTASK_

#define _MASSANDMOMENTUMTASK_

class BoundaryCondition;

#include "NairnMPM_Class/MPMTask.hpp"

class MassAndMomentumTask : public MPMTask
{
	public:
	
		// constructor
		MassAndMomentumTask(const char *);
	
		// required methods
		virtual void Execute(void);
	
	protected:
		CrackField cfld[2];
		double zDeriv[MaxShapeNds];
		
		void SetRigidBCs(long,int,double,double,BoundaryCondition **,BoundaryCondition **,BoundaryCondition **,BoundaryCondition **);
		void RemoveRigidBCs(BoundaryCondition **,BoundaryCondition **,BoundaryCondition **);
		void UnsetRigidBCs(BoundaryCondition **,BoundaryCondition **,BoundaryCondition **,BoundaryCondition **);
		void CountBCs(BoundaryCondition **,BoundaryCondition **,BoundaryCondition **);
};

#endif
