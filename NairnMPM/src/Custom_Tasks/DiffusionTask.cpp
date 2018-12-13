/********************************************************************************
    DiffusionTask.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Mon Mar 08 2004
    Copyright (c) 2003 John A. Nairn, All rights reserved.
 
    Diffusion calculations
   -------------------------
	See TransportTask.cpp comments with
		gTValue, gVCT, gQ in gDiff
 
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
#include "System/UnitsController.hpp"

// global
int DiffusionTask::active = NO_DIFFUSION;
double DiffusionTask::reference = 0.;				// zero-strain concentration
double DiffusionTask::viscosity = 0.001;			// poroelasticity viscosity
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
	
	// time step
	char fline[256];
	sprintf(fline,"time step (%s): %.7e",UnitsController::Label(ALTTIME_UNITS),transportTimeStep*UnitsController::Scaling(1.e3));
	if(active==POROELASTICITY_DIFFUSION)
		cout << "   Poroelasticity " << fline << endl;
	else
		cout << "   Diffusion " << fline << endl;
	cout << "   Time step factor: " << fmobj->GetTransCFLCondition() << endl;
	
	// featrues
	if(active==POROELASTICITY_DIFFUSION)
	{	cout << "   Reference pore pressure = " << reference*UnitsController::Scaling(1.e-6)
					<< " " << UnitsController::Label(PRESSURE_UNITS) << endl;
		cout << "   Fluid viscosity = " << viscosity*UnitsController::Scaling(1.e3)
					 << " " << UnitsController::Label(VISCOSITY_UNITS);
	}
	else
	{	char mline[100];
		sprintf(mline," =%8.4lf",reference);
		cout << "   Reference concentration" << mline;
	}
	cout << endl;
	
	return nextTask;
}

#pragma mark MASS AND MOMENTUM EXTRAPOLATIONS

// Task 1 Extrapolation of temperature to the grid
// Only called for non-rigid materials
TransportTask *DiffusionTask::Task1Extrapolation(NodalPoint *ndpt,MPMBase *mptr,double shape,short vfld,int matfld)
{
	double diffCT = theMaterials[mptr->MatID()]->GetDiffusionCT();
	double VpShape = mptr->GetVolume(DEFORMED_VOLUME)*diffCT*shape;
	double VpValueShape = mptr->pConcentration*VpShape;
	TransportField *gTrans = GetTransportFieldPtr(ndpt);
	gTrans->gTValue += VpValueShape;
	gTrans->gVCT += VpShape;
	Task1ContactExtrapolation(ndpt,vfld,matfld,VpValueShape,VpShape);
	return nextTask;
}

// Get Vp * CTp
double DiffusionTask::GetVpCTp(MPMBase *mptr)
{
	double diffCT = theMaterials[mptr->MatID()]->GetDiffusionCT();
	return mptr->GetVolume(DEFORMED_VOLUME)*diffCT;
}

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

// increment particle concentration with check in valid range
TransportTask *DiffusionTask::MoveTransportValue(MPMBase *mptr,double deltime,double rate) const
{
	mptr->pConcentration += deltime*rate;
	
	// limit concentration to 0 to 1
	if(mptr->pConcentration<0.)
		mptr->pConcentration = 0.;
	else if(mptr->pConcentration>1.)
		mptr->pConcentration = 1.;
	
	return nextTask;
}

#pragma mark ACCESSORS

// Return name of this task
const char *DiffusionTask::TaskName(void)
{	if(active==POROELASTICITY_DIFFUSION)
		return "pore pressure calculations";
	return "diffusion calculations";
}

// adjust time for given cell size if needed
TransportTask *DiffusionTask::TransportTimeStepFactor(int matid,double *diffCon)
{	*diffCon = theMaterials[matid]->MaximumDiffusion()/(theMaterials[matid]->GetDiffusionCT());
    return nextTask;
}

// return pointer on node to transport field
TransportField *DiffusionTask::GetTransportFieldPtr(NodalPoint *ndpt) const { return &(ndpt->gDiff); }

// return first boundary condition
NodalValueBC *DiffusionTask::GetFirstBCPtr(void) const { return firstConcBC; }
MatPtLoadBC *DiffusionTask::GetFirstFluxBCPtr(void) const { return firstFluxPt; }

// particle values
double *DiffusionTask::GetParticleValuePtr(MPMBase *mptr) const { return &(mptr->pConcentration); }
double *DiffusionTask::GetPrevParticleValuePtr(MPMBase *mptr) const { return &(mptr->pPreviousConcentration); }

// to check on diffusion or poroelasticity
bool DiffusionTask::HasDiffusion(void) { return active==MOISTURE_DIFFUSION; }
bool DiffusionTask::HasPoroelasticity(void) { return active==POROELASTICITY_DIFFUSION; }
bool DiffusionTask::HasFluidTransport(void) { return active!=NO_DIFFUSION; }

// convert poroelasticity MPa to Pa but no change to concentration potentions
double DiffusionTask::RescalePotential(void)
{	return active==POROELASTICITY_DIFFUSION ? UnitsController::Scaling(1.e6) : 1 ;
}

// convert Legacy poroelasticity (dV/V)/time to (dV/V)/time (factor=1) and
// convert concetration flux (kg/(m^2-sec) to Legacy (g/(mm^2 sec))
double DiffusionTask::RescaleFlux(void)
{	return active==POROELASTICITY_DIFFUSION ? 1. : UnitsController::Scaling(1.e-3);
}

