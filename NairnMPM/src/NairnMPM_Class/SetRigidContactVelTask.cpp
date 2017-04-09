/********************************************************************************
	SetRigidContactVelTask.cpp
	nairn-mpm-fea

	Created by John Nairn on May 12, 2016
	Copyright (c) 2016 John A. Nairn, All rights reserved.

	The tasks are:
	* Check all rigid contact materials (nmpmsNR to nmpmsRC)
	* Set the velocity only if a direction is controlled by a user-defined functions
	  (those not by function can only be constant velocity and no need to change)
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/SetRigidContactVelTask.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/RigidMaterial.hpp"

extern double mtime;

#pragma mark CONSTRUCTORS

SetRigidContactVelTask::SetRigidContactVelTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get mass matrix, find dimensionless particle locations,
//	and find grid momenta
void SetRigidContactVelTask::Execute(void)
{
    // Set rigid BC contact material velocities separately (so mass and momentum loop can be parallel)
	// GetVectorSetting() uses globals and therefore can't be parallel
	Vector newvel;
	bool hasDir[3];
//#pragma omp parallel for
	for(int p=nmpmsNR;p<nmpmsRC;p++)
	{   MPMBase *mpmptr = mpm[p];
		const RigidMaterial *matID = (RigidMaterial *)theMaterials[mpm[p]->MatID()];
		if(matID->GetVectorSetting(&newvel,hasDir,mtime,&mpmptr->pos))
		{   // change velocity if functions being used, otherwise keep velocity constant
			if(hasDir[0]) mpmptr->vel.x = newvel.x;
			if(hasDir[1]) mpmptr->vel.y = newvel.y;
			if(hasDir[2]) mpmptr->vel.z = newvel.z;
		}
	}
}
