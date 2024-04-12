/********************************************************************************
	PhaseFieldDiffusion.hpp
	nairn-mpm-fea
	
	Created by John Nairn on 2/17/2022
	Copyright (c) 2022 John A. Nairn, All rights reserved.
	
	Dependencies
		DiffusionTasks.hpp, TransportTask.hpp
 
    Phase field solved as a diffusion tasks. Tasks labeld "custom" need special code
   ----------------------------------------------------------------------------------
    Prepare for Calculations
        TransportTimeStepFactor() - for for time step (custom to call MaximumPhaseDiffusivity())
        Initialize() - print details (custom)
    Initialization:
        InitializeForTimeStep() - built-in zeros all diffusion gTValue, gVCT, and gQ on node
    Mass and Momentum Task
        Task1Extrapolation() - extrapolate gTValue, gVCT using viscosity (custom)
        Task1Reduction() - uadd ghost to real
    Post M&M Extrapolation Task
        GetTransportNodalValue() - divide gTValue by gVCT
        TransportBCsAndGradients() - custom is GetFirstBCPtr() for this task only
            ImposeValueBCs(): - impose grid BCs
            GetGradients(): standard methods (based on diffusion task number)
            RestoreValueBCs() - paste back initial gTValue (in FMPM version always)
    Update strains first (if used)
        NEEDED: to support other then USL, need to scale dPhase (or do all at one time?)
        Check on transport properties usages
    Grid Forces Task
        AddForces():Add transport force and phase source terms to gQ (custom)
        ForcesReduction(), add gQ from ghost to real
    Post Forces Task
        SetTransportFluxBCs() - impose flux BCs in gQ (custom is GetFirstFluxBCPtr() for this task only)
            (In transport only, this is done in GridForcesTask.cpp and PostForcesTask.cpp is not used)
    Update Momenta Task
        UpdateTransport() - divide gQ by gVCT to get rates and update nodal value
        ImposeValueGridBCs() - impose grid BCs (custom is GetFirstBCPtr() for this task only)
    NEEDED (check) - XPIC/FMPM Calculations
        Calculate gTstar = Sum gTk starting with gTk = m*gTValue   (* to fix comment formatting)
    Update Particles Task
        If FMPM order>1: Copy gTstar to gTValue (if order>1 in FMPM)
            Impose grid BCs (ImposeValueBCs(time,false))
        Particle loop:
            Extrapolate value and rate to particle (IncrementValueExtrap() and IncrementTransportRate())
            AdjustRateAndValue() - keep 0 to 1 and rate positive (custom, phase field materials only)
            GetDeltaValue() - get phase change based on grid extrapolations
                (custom, phase field materials only and save in history)
            MoveTransportValue() - update particle value (custom - phase file only and save on history)
    Update strains last (if used and if second extrapolation)
        NEEDED: to support other then USL, need to scale dPhase (or do all at one time?)
        Check on transport properties usages

    Material Point class support
        pDiff[number] - stores DiffusionInfo structure (conc, prevConc, grad)
    Nodal Point class support
        gDiff[number] -  stores TransportField (for extrapolation and XPIC stuff)
    Material class support
        SupportsPhaseField() - true or false if phase field for this task; rest only called after this known to be true
            MaximumPhaseDiffusivity() - for time step
            GetPhaseFieldViscosity() - viscosity
            GetPhaseFieldDTensor() -  "diffusion" tensor tensor
			GetMatDiffusionSource() - return phase-field defined source term (use this and/or the next for source terms)
 			GetParticleDiffusionSource() - return phase-field defined source term depending on particle state
            AdjustPhaseFieldValue() - adjust extrapolated value and rate if needed
            StorePhaseFieldValue() - store phase field value
            StorePhaseFieldDelta() - store phase field delta value
        By default, particle deletion sets all pDiff[i] of other diffusion and all history doubles to zero
            A material that needs something different can override ResetHistoryData() to set history and pDiff[i]
    Boundary conditions
        Main code divides nodal and flux BCs into lists starting at firstDiffBC[i] and firstDiffFluxBC[i]
        Diffusion tasks returns the right BC depending on phase style parameters
********************************************************************************/

#include "stdafx.h"
#include "Custom_Tasks/PhaseFieldDiffusion.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "System/UnitsController.hpp"

#pragma mark CONSTRUCTORS

// Constructors
// diffStyle can support various types of phase field diffusion at the same time)
PhaseFieldDiffusion::PhaseFieldDiffusion(int diffStyle) : DiffusionTask(0.,1.,diffStyle)
{
}

#pragma mark MASS AND MOMENTUM EXTRAPOLATIONS

// Task 1 Extrapolation of phase field to the grid
TransportTask *PhaseFieldDiffusion::Task1Extrapolation(NodalPoint *ndpt,MPMBase *mptr,double shape,short vfld,int matfld)
{
    // skip if not supported
    MaterialBase *matref = theMaterials[mptr->MatID()];
    if(!matref->SupportsPhaseField(style)) return nextTask;

	// gTValue and gVCT
    double viscosity = matref->GetPhaseFieldViscosity(style);
	double VpShape = mptr->GetVolume(DEFORMED_VOLUME)*shape;
	double VpEtaShape = VpShape*viscosity;
	double pPhase = GetParticleValue(mptr);
	TransportField *gTrans = GetTransportFieldPtr(ndpt);
	gTrans->gTValue += pPhase*VpEtaShape;
	gTrans->gVCT += VpEtaShape;

	return nextTask;
}

// Get Vp * CTp (used in XPIC extrapolations)
double PhaseFieldDiffusion::GetVpCTp(MPMBase *mptr)
{   // skip if not supported
    MaterialBase *matref = theMaterials[mptr->MatID()];
    if(!matref->SupportsPhaseField(style)) return 0.;
    
    // Vp * eta
	return mptr->GetVolume(DEFORMED_VOLUME)*matref->GetPhaseFieldViscosity(style);
}

#pragma mark GRID FORCES EXTRAPOLATIONS

// find forces for phase field calculation (non-rigid particles only)
TransportTask *PhaseFieldDiffusion::AddForces(NodalPoint *ndptr,MPMBase *mptr,double sh,double dshdx,
                        double dshdy,double dshdz,TransportProperties *t,short vfld,int matfld)
{
	// change to phase field diffusion tensor (skip if not one)
    MaterialBase *matref = theMaterials[mptr->MatID()];
    if(!matref->SupportsPhaseField(style)) return nextTask;
    
    // get custom diffusion tensor
	TransportProperties tPhase;
	tPhase.diffusionTensor = matref->GetPhaseFieldDTensor(style,mptr);
	
	// internal force using custom diffusion tensor
    TransportField *gTrans = GetTransportFieldPtr(ndptr);
    gTrans->gQ += mptr->FDiff(dshdx,dshdy,dshdz,&tPhase,number);
	
    // Phase field source terms
	gTrans->gQ += matref->GetMatDiffusionSource(style,mptr,GetParticleValue(mptr),
								mptr->GetVolume(DEFORMED_VOLUME),sh,dshdx,dshdy,dshdz);

	return nextTask;
}

#pragma mark UPDATE PARTICLES TASK

// after extrapolated, give task a chance to change rate and value if needed
void PhaseFieldDiffusion::AdjustRateAndValue(MPMBase *mptr,double &value,
											 double &rate,double &lumpedValue,double deltime) const
{
	// skip if not phase field material
    MaterialBase *matref = theMaterials[mptr->MatID()];
	if(!matref->SupportsPhaseField(style)) return;
	
    // save on particle too (if changed, changed lumpdValue too)
    matref->AdjustPhaseFieldValue((DiffusionTask *)this,GetParticleValue(mptr),
						GetPrevParticleValue(mptr),value,rate,lumpedValue,deltime);
}

// increment particle concentration (time is always timestep)
TransportTask *PhaseFieldDiffusion::MoveTransportValue(MPMBase *mptr,double deltime,double rate,double value) const
{
    // skip if not phase field material
    MaterialBase *matref = theMaterials[mptr->MatID()];
    if(!matref->SupportsPhaseField(style)) return nextTask;
	
	// update particle value (in pDiff[number])
	TransportTask::MoveTransportValue(mptr,deltime,rate,value);
	
	// save on particle too
    matref->StorePhaseFieldValue(style,mptr,GetParticleValue(mptr));

	return nextTask;
}

#pragma mark UPDATE PARTICLE STRAIN TASK

// after extrapolated, find change this update on particle extrapolated from the grid
//  Do not set *dV, get local change and store in history variable
void PhaseFieldDiffusion::GetDeltaValue(MPMBase *mptr,double pValueExtrap,double *dV) const
{	// skip if not phase field materials
	MaterialBase *matref = theMaterials[mptr->MatID()];
	if(!matref->SupportsPhaseField(style)) return;
	
	// get change in phase field from extrapolated values
	double dPhase;
	TransportTask::GetDeltaValue(mptr,pValueExtrap,&dPhase);
	
    // save delta on particle here and not passed back in *dV
    matref->StorePhaseFieldDelta(style,mptr,dPhase);
}

#pragma mark ACCESSORS

// Return style of this diffusion this task
const char *PhaseFieldDiffusion::StyleName(void)
{	switch(style)
	{   case FRACTURE_PHASE_FIELD:
		   return "fracture phase field";
		   break;
	   case BATTERY_PHASE_FIELD:
		   return "battery phase field";
		   break;
	   case CONDUCTION_PHASE_FIELD:
		   return "electrical conduction";
		   break;
	   default:
		   break;
   }
	return DiffusionTask::StyleName();
}

// adjust time for given cell size if needed
TransportTask *PhaseFieldDiffusion::TransportTimeStepFactor(int matid,double *diffCon)
{   if(theMaterials[matid]->SupportsPhaseField(style))
        *diffCon = theMaterials[matid]->MaximumPhaseDiffusivity(style);
    else
        *diffCon = 0.;
    return nextTask;
}

