/********************************************************************************
	UpdateStrainsFirstTask.cpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	Update strains on particles after initial projection to the grid
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
	FullStrainUpdate(strainTimestep,FALSE,fmobj->np);
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
	
	// loop over nonrigid particles
	// This works as parallel when material properties change with particle state because
	//	all such materials should create a copy of material properties in the threads
//#pragma omp parallel for
	for(int p=0;p<nmpmsNR;p++)
    {   // next particle
        MPMBase *mptr = mpm[p];
        
        // this particle's material
        const MaterialBase *matRef = theMaterials[mptr->MatID()];
        
        // make sure have mechanical properties for this material and angle
		// Must replace with get copy of mechanical properties
		void *properties = matRef->GetCopyOfMechanicalProps(mptr,np);
        
        // finish on the particle
		mptr->UpdateStrain(strainTime,secondPass,np,properties,matRef->GetField());
		
		// delete properties
		matRef->DeleteCopyOfMechanicalProps(properties,np);
    }
}
