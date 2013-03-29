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
#include "Nodes/NodalPoint.hpp"

#include "Materials/MaterialBase.hpp"

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

	FullStrainUpdate(strainTimestep,FALSE,fmobj->np);
	
#ifdef _PROFILE_TASKS_
	totalTaskTime+=fmobj->CPUTime()-beginTime;
#endif
}

#pragma mark MPMBase Class Methods

/**********************************************************
 Update Strains on all particles
 Must impose grid velocity BCs and velocity
 alterations due to contact first
 secondPass will be TRUE only for USAVG method
 **********************************************************/
void UpdateStrainsFirstTask::FullStrainUpdate(double strainTime,int secondPass,int np)
{
    NodalPoint::GetGridVelocitiesForStrainUpdate();			// velocities needed for strain update
	
	// loop over non rigid particles
	for(MPMBase::currentParticleNum=0;MPMBase::currentParticleNum<nmpmsNR;MPMBase::currentParticleNum++)
    {   // next particle
        MPMBase *mptr = mpm[MPMBase::currentParticleNum];
        
        // this particle's material
        MaterialBase *matRef=theMaterials[mptr->MatID()];
        
        // exit if rigid (in case some before the last non rigid one)
        if(matRef->Rigid()) return;
        
        // make sure have mechanical properties for this material and angle
        matRef->LoadMechanicalProps(mptr,np);
        
        // finish on the particle
		mptr->UpdateStrain(strainTime,secondPass,np,matRef,matRef->GetField());
    }
}
