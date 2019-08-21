/********************************************************************************
	ProjectRigidBCsTask.hpp
	nairn-mpm-fea

	Created by John Nairn on March 6, 2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _PROJECTRIGIDBCSTASK_

#define _PROJECTRIGIDBCSTASK_

class BoundaryCondition;

#include "NairnMPM_Class/MPMTask.hpp"

class ProjectRigidBCsTask : public MPMTask
{
	public:
	
		// constructor
		ProjectRigidBCsTask(const char *);
	
		// required methods
		virtual void Execute(int);
	
		// class methods
		static void UnsetRigidBCs(BoundaryCondition **,BoundaryCondition **,BoundaryCondition **,BoundaryCondition **);
		static void SetRigidBCs(int,int,int,double,double,int,BoundaryCondition **,BoundaryCondition **,BoundaryCondition **,BoundaryCondition **);
		static void RemoveRigidBCs(BoundaryCondition **,BoundaryCondition **,BoundaryCondition **);
		static void CountBCs(BoundaryCondition **,BoundaryCondition **,BoundaryCondition **);

	protected:
};

#endif
