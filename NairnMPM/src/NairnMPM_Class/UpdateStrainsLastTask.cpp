/********************************************************************************
	UpdateStrainsLastTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Update strains on particles after they have been updated. This task
	is active for USL and USAVG update methods
 
	The tasks are:
	--------------
    * Update nodal transport values
	* Full strain update
		- Get grid velocities (p/m)
		- Copy mechanical properties of a material
		- Tell material to update strain (etc) on the particle
		- Also extrapolates transport to particle and saves as "previous"
 
	The tasks for transport only are:
	---------------------------------
	* Update nodal transport values
	* Full strain update
		- Tell material to particle heat energy
		- Also extrapolates transport to particle and saves as "previous"
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/UpdateStrainsLastTask.hpp"
#include "NairnMPM_Class/UpdateStrainsFirstTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Custom_Tasks/TransportTask.hpp"

#pragma mark CONSTRUCTORS

UpdateStrainsLastTask::UpdateStrainsLastTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get total grid point forces (except external forces)
void UpdateStrainsLastTask::Execute(int taskOption)
{
	// update strains based on current velocities
	UpdateStrainsFirstTask::FullStrainUpdate(strainTimestepLast,(fmobj->mpmApproach==USAVG_METHOD),fmobj->np,true);
}
