/********************************************************************************
	UpdateStrainsFirstTask.cpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	Update strains on particles after initial projection to the grid
 
	Input Variables
		mpm[]->ncpos
		nd[]->vel, gTemperature, gConcentration
 
	Output Variables
		theMaterials[]->LoadMechanicalProps() - changes any properties that depend
												on particle state
		MPMBase::currentParticleNum - used in strain update loop
		mpm[]->sp, ep, eplast, wrot, plastEnergy, dispEnergy, strainEnergy,
				extWork, matData->, pPreviousTemperature, pPreviousConcentration
		ConductionTask::dTemperature
		DiffusionTask::dConcentration

********************************************************************************/

#include "NairnMPM_Class/UpdateStrainsFirstTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "MPM_Classes/MPMBase.hpp"

#pragma mark CONSTRUCTORS

UpdateStrainsFirstTask::UpdateStrainsFirstTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Update strains with just-extrapolated momenta
void UpdateStrainsFirstTask::Execute(void)
{
#ifdef _PROFILE_TASKS_
	double beginTime=fmobj->CPUTime();
#endif

	MPMBase::FullStrainUpdate(strainTimestep,FALSE,fmobj->np);
	
#ifdef _PROFILE_TASKS_
	totalTaskTime+=fmobj->CPUTime()-beginTime;
#endif
}
