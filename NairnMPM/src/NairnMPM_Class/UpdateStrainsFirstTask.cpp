/********************************************************************************
	UpdateStrainsFirstTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	Update strains on particles after initial projection to the grid
 
	The tasks are:
	--------------
	* Full strain update
		- Get grid velocities (p/m)
		- Copy mechanical properties of a material
		- Tell material to update strain (etc) on the particle
		- Also extrapolates transport to particle and saves a "previous"

	The tasks for transport only are:
	---------------------------------
	* Full strain update
		- Tell material to particle heat energy
		- Also extrapolates transport to particle and saves a "previous"
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/UpdateStrainsFirstTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Materials/MaterialBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "NairnMPM_Class/XPICExtrapolationTask.hpp"
#include "Global_Quantities/BodyForce.hpp"

// class globals
UpdateStrainsFirstTask *USFTask = NULL;

// one copy of property buffers for each thread
void **UpdateStrainsFirstTask::matBuffer = NULL;
void **UpdateStrainsFirstTask::altBuffer = NULL;

#pragma mark CONSTRUCTORS

UpdateStrainsFirstTask::UpdateStrainsFirstTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Update strains with just-extrapolated momenta
bool UpdateStrainsFirstTask::Execute(int taskOption)
{
	FullStrainUpdate(strainTimestepFirst,false,fmobj->np,false);
    return true;
}

#ifdef RESTART_OPTION
bool UpdateStrainsFirstTask::BlockRestartOption(void) const { return true; }
#endif

#pragma mark UpdateStrainFirstTask Class Methods

// create buffers large enough for all active materials and plastic laws
// need one for each thread
// throws CommonException()
void UpdateStrainsFirstTask::CreatePropertyBuffers(int numThreads)
{
    matBuffer = new (nothrow) void *[numThreads];
    if(matBuffer==NULL)
        throw CommonException("Out of memory allocating memory for propery buffers","UpdateStrainsFirstTask::CreatePropertyBuffers");
    altBuffer = new (nothrow) void *[numThreads];
    if(altBuffer==NULL)
        throw CommonException("Out of memory allocating memory for propery buffers","UpdateStrainsFirstTask::CreatePropertyBuffers");
    
    for(int i=0;i<numThreads;i++)
    {   if(MaterialBase::maxPropertyBufferSize>0)
        {   matBuffer[i] = (void *)(new (nothrow) char[MaterialBase::maxPropertyBufferSize]);
            if(matBuffer[i]==NULL)
                throw CommonException("Out of memory allocating memory for propery buffers","UpdateStrainsFirstTask::CreatePropertyBuffers");
        }
        else
            matBuffer[i]=NULL;
        if(MaterialBase::maxAltBufferSize>0)
        {   altBuffer[i] = (void *)(new(nothrow) char[MaterialBase::maxAltBufferSize]);
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
	throws CommonException()
**********************************************************/
void UpdateStrainsFirstTask::FullStrainUpdate(double strainTime,int secondPass,int np,bool postUpdate)
{
	CommonException *usfErr = NULL;

	// Get velocities on the grid for extrapolating velocity gradient
	int xpicUsed = bodyFrc.UsingFMPM() ? bodyFrc.GetXPICOrder() : 0;
	if(xpicUsed>1)
	{	// FMPM(k>1) only here - do XPIC calculations if needed
		if(!postUpdate || !fmobj->skipPostExtrapolation)
		{	// !postUpdate means USF, USAVG-, or USAVG+ - always do calcualtions
			// postUpdate with !skipPostExtrapolation  means USAVG+ or USL+ - needs new calculations
			// postUpdate with skipPostExtrapolation  means USAVG- or USL- - use stored results
			XPICMechanicsTask->Execute(0);
			// velocity BCs
			if(bodyFrc.GridBCOption()!=GRIDBC_LUMPED_ONLY)
				NodalVelBC::GridVelocityConditions(UPDATE_GRID_STRAINS_CALL);
		}
	}
	else
	{	// FLIP, XPIC(k), and FMPM(1) use lumped mass matrix
#pragma omp parallel for
		for(int i=1;i<=*nda;i++)
			nd[nda[i]]->GridValueCalculation(VELOCITY_FOR_STRAIN_UPDATE);
	}

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
		catch(CommonException& err)
		{	if(usfErr==NULL)
			{
#pragma omp critical (error)
				usfErr = new CommonException(err);
			}
		}
		catch(std::bad_alloc&)
		{	if(usfErr==NULL)
			{
#pragma omp critical (error)
				usfErr = new CommonException("Memory error","UpdateStrainsFirstTask::FullStrainUpdat");
			}
		}
		catch(...)
		{	if(usfErr==NULL)
			{
#pragma omp critical (error)
				usfErr = new CommonException("Unexpected error","UpdateStrainsFirstTask::FullStrainUpdat");
			}
		}
    }
	
	// throw error if it occurred
	if(usfErr!=NULL) throw *usfErr;
}
