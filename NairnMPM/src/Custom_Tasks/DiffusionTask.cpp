/********************************************************************************
    DiffusionTask.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Mon Mar 08 2004
    Copyright (c) 2003 John A. Nairn, All rights reserved.
 
    Diffusion calculations
   -------------------------
	See TransportTask.cpp comments with
		gTValue, gMTp, gQ in gDiff
 
    Update Particles Task
        cut off particle potential to range 0 to 1
    Chemical potential (or concentration potential 0 to 1)
        Internally all calculations in terms or potential (0 to 1)
        Output concentration and concentration gradient scaled to
            csat to get weight fraction and weight fraction gradient
        Flux set as mass per area per sec. See documentation for conversions.
********************************************************************************/

#include "stdafx.h"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "Boundary_Conditions/NodalConcBC.hpp"
#include "Boundary_Conditions/MatPtFluxBC.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Exceptions/CommonException.hpp"

// global
bool DiffusionTask::active=FALSE;
double DiffusionTask::reference = 0.;				// zero-strain concentration
DiffusionTask *diffusion=NULL;

#pragma mark STANDARD METHODS

// called once at start of MPM analysis and after preliminary calcse are eon
TransportTask *DiffusionTask::Initialize(void)
{
	// allocate diffusion data on each particle
    // done before know number of nonrigid, so do on all
	for(int p=0;p<nmpms;p++)
	{	mpm[p]->pDiffusion = new double[3];
		for(int i=0;i<3;i++) mpm[p]->pDiffusion[i] = 0.;
	}
	
	cout << "Coupled " << TaskName() << endl;
	char mline[100];
	sprintf(mline,"   Reference concentration =%8.4lf",reference);
	cout << mline << endl;
	
	return nextTask;
}

#pragma mark MASS AND MOMENTUM EXTRAPOLATIONS

// zero gradients on the particle
void DiffusionTask::ZeroTransportGradients(MPMBase *mptr)
{	mptr->AddConcentrationGradient();
}

// Add gradients on the particles
void DiffusionTask::AddTransportGradients(MPMBase *mptr,Vector *deriv,NodalPoint *ndptr,short vfld)
{	mptr->AddConcentrationGradient(ScaleVector(deriv,ndptr->gDiff.gTValue));
}

#pragma mark GRID FORCES EXTRAPOLATIONS

// find forces for diffusion calculation (mm^3/sec) (non-rigid particles only)
TransportTask *DiffusionTask::AddForces(NodalPoint *ndptr,MPMBase *mptr,double sh,double dshdx,
										double dshdy,double dshdz,TransportProperties *t,short vfld,int matfld)
{	// internal force
	ndptr->gDiff.gQ += mptr->FDiff(dshdx,dshdy,dshdz,t);
	
	// add source terms (should be potential per sec, if c units per second, divide by concSaturation)
	
	return nextTask;
}

#pragma mark UPDATE PARTICLES TASK

// increment particle concentration with check in valud range
TransportTask *DiffusionTask::MoveTransportValue(MPMBase *mptr,double deltime,double rate) const
{	mptr->pConcentration += deltime*rate;
    if(mptr->pConcentration<0.)
        mptr->pConcentration = 0.;
    else if(mptr->pConcentration>1.)
        mptr->pConcentration = 1.;
	return nextTask;
}

#pragma mark ACCESSORS

// Return name of this task
const char *DiffusionTask::TaskName(void) { return "diffusion calculations"; }

// adjust time for given cell size if needed
TransportTask *DiffusionTask::TransportTimeStepFactor(int matid,double *diffCon)
{	*diffCon = theMaterials[matid]->MaximumDiffusion();
    return nextTask;
}

// Get transport mass, mTp in notes, which here is V_p and get particle value
double DiffusionTask::GetTransportMassAndValue(MPMBase *mptr,double *pTValue)
{	*pTValue = mptr->pConcentration;
    return mptr->GetVolume(DEFORMED_VOLUME);
}

// return point on node to transport field
TransportField *DiffusionTask::GetTransportFieldPtr(NodalPoint *ndpt) const { return &(ndpt->gDiff); }

// return first boundary condition
NodalValueBC *DiffusionTask::GetFirstBCPtr(void) const { return firstConcBC; }
MatPtLoadBC *DiffusionTask::GetFirstFluxBCPtr(void) const { return firstFluxPt; }

// particle values
double *DiffusionTask::GetParticleValuePtr(MPMBase *mptr) const { return &(mptr->pConcentration); }
double *DiffusionTask::GetPrevParticleValuePtr(MPMBase *mptr) const { return &(mptr->pPreviousConcentration); }

		
