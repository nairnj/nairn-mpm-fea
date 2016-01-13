/********************************************************************************
	UpdateStrainsFirstTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	Update strains on particles after initial projection to the grid
********************************************************************************/

#include "NairnMPM_Class/UpdateStrainsFirstTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Materials/MaterialBase.hpp"
#include "Exceptions/CommonException.hpp"

// one copy of property buffers for each thread
void **UpdateStrainsFirstTask::matBuffer = NULL;
void **UpdateStrainsFirstTask::altBuffer = NULL;

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

#pragma mark UpdateStrainFirstTask Class Methods

// create buffers large enough for all active materials and plastic laws
// need one for each thread
void UpdateStrainsFirstTask::CreatePropertyBuffers(int numThreads)
{
    matBuffer = (void **)malloc(numThreads*sizeof(void *));
    if(matBuffer==NULL)
        throw CommonException("Out of memory allocating memory for propery buffers","UpdateStrainsFirstTask::CreatePropertyBuffers");
    altBuffer = (void **)malloc(numThreads*sizeof(void *));
    if(altBuffer==NULL)
        throw CommonException("Out of memory allocating memory for propery buffers","UpdateStrainsFirstTask::CreatePropertyBuffers");
    
    for(int i=0;i<numThreads;i++)
    {   if(MaterialBase::maxPropertyBufferSize>0)
        {   matBuffer[i] = (void *)malloc(MaterialBase::maxPropertyBufferSize);
            if(matBuffer[i]==NULL)
                throw CommonException("Out of memory allocating memory for propery buffers","UpdateStrainsFirstTask::CreatePropertyBuffers");
        }
        else
            matBuffer[i]=NULL;
        if(MaterialBase::maxAltBufferSize>0)
        {   altBuffer[i] = (void *)malloc(MaterialBase::maxAltBufferSize);
            if(altBuffer[i]==NULL)
                throw CommonException("Out of memory allocating memory for propery buffers","UpdateStrainsFirstTask::CreatePropertyBuffers");
        }
        else
            altBuffer[i]=NULL;
    }
}

/**********************************************************
    Update Strains on all particles
    Must impose grid velocity BCs and velocity
        alterations due to contact first
    secondPass will be TRUE only for USAVG method
**********************************************************/
void UpdateStrainsFirstTask::FullStrainUpdate(double strainTime,int secondPass,int np)
{
	CommonException *usfErr = NULL;
	
    NodalPoint::GetGridVelocitiesForStrainUpdate();			// velocities needed for strain update

	// loop over nonrigid particles
	// This works as parallel when material properties change with particle state because
	//	all such materials should create a copy of material properties in the threads
#pragma omp parallel for
	for(int p=0;p<nmpmsNR;p++)
    {   // next particle
        MPMBase *mptr = mpm[p];
        
        // this particle's material
        const MaterialBase *matRef = theMaterials[mptr->MatID()];
        
		try
		{	// make sure have mechanical properties for this material and angle
            int tn = GetPatchNumber();
			void *properties = matRef->GetCopyOfMechanicalProps(mptr,np,matBuffer[tn],altBuffer[tn],0);
			
			// finish on the particle
			mptr->UpdateStrain(strainTime,secondPass,np,properties,matRef->GetField());
		}
		catch(CommonException err)
		{	if(usfErr==NULL)
			{
#pragma omp critical (error)
				usfErr = new CommonException(err);
			}
		}
    }
	
	// throw error if it occurred
	if(usfErr!=NULL) throw *usfErr;
}
