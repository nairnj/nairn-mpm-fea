/********************************************************************************
    DiffusionTask.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Mon Mar 08 2004
    Copyright (c) 2003 John A. Nairn, All rights reserved.
 
    Diffusion calculations
   -------------------------
	See TransportTask.cpp comments with
		gTValue, gVCT, gQ in gDiff[]
		(materials with csat depending on position extrapolate gTValueRel,
			or gTRelValue and uses one of them in AddTransportGradients.
 
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
#ifdef TRANSPORT_FMPM
    #include "NairnMPM_Class/XPICExtrapolationTaskTO.hpp"
#endif

// globals
DiffusionTask *diffusion=NULL;
DiffusionTask *otherDiffusion=NULL;
int numDiffusion=-1;

#pragma mark CONSTRUCTORS

// Constructors
DiffusionTask::DiffusionTask(double prop1,double prop2,int diffusionStyle) : TransportTask()
{
	reference = prop1;
	viscosity = prop2;			// only for poroelasticity
	style = diffusionStyle;     // 1 for diffusion, 2 forporoelasticity, id for subclass PhaseFieldDiffusion
	number = 0;					// for diffusion, others will get a number when CountDiffusionTasks() called
}

#pragma mark STANDARD METHODS

// called once at start of MPM analysis and after preliminary calcs are done
TransportTask *DiffusionTask::Initialize(void)
{
	cout << "Coupled " << TaskName() << " (number " << number << ")" << endl;
	cout << "   For " << StyleName() << " (ID " << style << ")" << endl;
	
	// time step
	char fline[256];
	size_t fsize=256;
	snprintf(fline,fsize,"   Time step maximum (%s): %.7e",UnitsController::Label(ALTTIME_UNITS),
								transportTimeStep*UnitsController::Scaling(1.e3));
	cout << fline << endl;
	cout << "   Time step factor: " << fmobj->GetTransCFLCondition() << endl;
	
	// features
	if(style==POROELASTICITY_DIFFUSION)
	{	cout << "   Reference pore pressure = " << reference*UnitsController::Scaling(1.e-6)
					<< " " << UnitsController::Label(PRESSURE_UNITS) << endl;
		cout << "   Fluid viscosity = " << viscosity*UnitsController::Scaling(1.e3)
					 << " " << UnitsController::Label(VISCOSITY_UNITS) << endl;
	}
	else if(style==MOISTURE_DIFFUSION)
	{	snprintf(fline,fsize," =%8.4lf",reference);
		cout << "   Reference concentration" << fline << endl;
	}
	
	return nextTask;
}

#pragma mark MASS AND MOMENTUM EXTRAPOLATIONS

// Task 1 Extrapolation of concentration or pore pressure to the grid
// Only called for non-rigid materials
TransportTask *DiffusionTask::Task1Extrapolation(NodalPoint *ndpt,MPMBase *mptr,double shape,short vfld,int matfld)
{
	MaterialBase *matref = theMaterials[mptr->MatID()];
	double diffCT = matref->GetDiffusionCT();
	double pConc = GetParticleValue(mptr);
	double VpCTShape = mptr->GetVolume(DEFORMED_VOLUME)*diffCT*shape;
	double VpValueShape = pConc*VpCTShape;
	
	TransportField *gTrans = GetTransportFieldPtr(ndpt);
	gTrans->gTValue += VpValueShape;
	gTrans->gVCT += VpCTShape;
	
	if(doCrelExtrapolation)
	{	// only true for MOISTURE_DIFFUSION when some material varies Csat with position
		//   and it allows saturation to change with particle state
		// if csatRelative<1, pConc may exceed 1
		double csatRel = matref->GetCsatRelative(mptr);
#ifdef USE_GTVALUEREL
		gTrans->gTValueRel += VpValueShape/csatRel;
#else
		gTrans->gTRelValue += VpCTShape*csatRel;
#endif
	}

	
	Task1ContactExtrapolation(ndpt,vfld,matfld,VpValueShape,VpCTShape);
	return nextTask;
}

// Get Vp * CTp
double DiffusionTask::GetVpCTp(MPMBase *mptr)
{   double diffCT = theMaterials[mptr->MatID()]->GetDiffusionCT();
	return mptr->GetVolume(DEFORMED_VOLUME)*diffCT;
}

// zero gradients on the particle
void DiffusionTask::ZeroTransportGradients(MPMBase *mptr)
{	mptr->AddConcentrationGradient(number);
}

// Add gradients on the particles
void DiffusionTask::AddTransportGradients(MPMBase *mptr,Vector *deriv,NodalPoint *ndptr,short vfld)
{
#ifdef USE_GTVALUEREL
	double theValue = !doCrelExtrapolation ? ndptr->gDiff[number].gTValue :
												ndptr->gDiff[number].gTValueRel ;
#else
	double theValue = !doCrelExtrapolation ? ndptr->gDiff[number].gTValue :
						ndptr->gDiff[number].gTValue/ndptr->gDiff[number].gTRelValue;
#endif
	mptr->AddConcentrationGradient(number,ScaleVector(deriv,theValue));
}

#pragma mark GRID FORCES EXTRAPOLATIONS

// find forces for diffusion calculation (mm^3/sec) (non-rigid particles only)
TransportTask *DiffusionTask::AddForces(NodalPoint *ndptr,MPMBase *mptr,double sh,double dshdx,
										double dshdy,double dshdz,TransportProperties *t,short vfld,int matfld)
{
	// the material point's material
	MaterialBase *matref = theMaterials[mptr->MatID()];
	
	// internal force
	if(!doCrelExtrapolation)
		ndptr->gDiff[number].gQ += mptr->FDiff(dshdx,dshdy,dshdz,t,number);
	else
	{	// scale by D by csatRalative on the particle (never for poroelasticity)
		double scale = matref->GetCsatRelative(mptr);
		TransportProperties tPhase = *t;
		ScaleTensor(&tPhase.diffusionTensor,scale);
		ndptr->gDiff[number].gQ += mptr->FDiff(dshdx,dshdy,dshdz,&tPhase,number);
	}
	
	// All materials to add a source term
	// Normally the return is Vp*sh*source, but gradient passed in case needed
	ndptr->gDiff[number].gQ += matref->GetMatDiffusionSource(style,mptr,GetParticleValue(mptr),
								mptr->GetVolume(DEFORMED_VOLUME),sh,dshdx,dshdy,dshdz);

	return nextTask;
}

#pragma mark UPDATE PARTICLES TASK

// after extrapolated, give task a chance to change rate and value(s) if needed
// only called for moisture and poroelasticity
void DiffusionTask::AdjustRateAndValue(MPMBase *mptr,double &value,
									   double &rate,double &lumpedValue,double deltime) const
{
	// Poroelasticity does not change anything
	if(style==POROELASTICITY_DIFFUSION) return;
	
#ifdef NO_UPPER_LIMIT
	// note value must be above zero
	if(value<0.) value = 0.;
	
	// do lumped value just in case (only needed for Blended FLIP/FMPM(k>1)
	if(lumpedValue<0.) lumpedValue = 0.;
	
	// check rate, unless pure FMPM
	if(!usingXPIC || usingFraction<1.)
	{	// grab particle value
		double pConc = GetParticleValue(mptr);

		// make sure rate does not jump negative (only needed for FLIP)
		// minimum rate is when pConc + mincdt*deltime = 0 or mincdt = -pConc/deltim
		double mincdt = -pConc/deltime;
		if(rate<mincdt) rate = mincdt;
	}
#else
	// new value must stay between 0 and csatRelative for both using XPIC or FLIP
	double csatRelative = 1.;
	if(doCrelExtrapolation)
	{	MaterialBase *matref = theMaterials[mptr->MatID()];
		csatRelative = matref->GetCsatRelative(mptr);
	}
	value = fmax(fmin(csatRelative,value),0.);
	
	// do lumped value just in case (only needed for Blended FLIP/FMPM(k>1)
	lumpedValue = fmax(fmin(csatRelative,lumpedValue),0.);
	
	// check rate, unless pure FMPM
	if(!usingXPIC || usingFraction<1.)
	{	// grab particle value
		double pConc = GetParticleValue(mptr);
		
		// make sure rate does not jump outside the range
		// minimum rate for FLIP is when pConc+mincdt*dt=0 or mincdt = -pConc/dt
		double mincdt = -pConc/deltime;
		if(rate<mincdt)
		{	rate = mincdt;
		}
		else
		{	// maximum rate is when pConc+maxcdt*dt = csatRelative
			double maxcdt = csatRelative/deltime+mincdt;
			if(rate>maxcdt)
			{	rate = maxcdt;
			}
		}
	}
#endif
}

// After extrapolate, find dC value in super class method
// Pororelasticity prevent prevConc from going below zere if if needed
//		sets dC to match.
void DiffusionTask::GetDeltaValue(MPMBase *mptr,double pValueExtrap,double *dC) const
{	// get change in phase field from extrapolated values
	TransportTask::GetDeltaValue(mptr,pValueExtrap,dC);
	
	// Poroelasticity prevents pPrev0 from going below zero
	if(style==POROELASTICITY_DIFFUSION)
	{	if(mptr->pDiff[0]->prevConc<0.)
		{	// if extrpolated value+dpud has gone negative, then bring
			// it back to zero and change dC to reach zero
			// Note: mptr->pDiff[0]->prevConc-res.dC = pPrevConc0
			*dC = -fmax(mptr->pDiff[0]->prevConc-*dC,0.);
			mptr->pDiff[0]->prevConc = 0.;
		}
	}
}

// increment particle concentration with check in valid range
TransportTask *DiffusionTask::MoveTransportValue(MPMBase *mptr,double deltime,double rate,double value) const
{
	// grab pointer to value
	double *pConc = GetParticleValuePtr(mptr);
	
#ifdef TRANSPORT_FMPM
	if(usingXPIC)
	{	// FMPM update just replace the value or blend with FLIP
		if(usingFraction<1.)
			*pConc = (1.-usingFraction)*(*pConc+deltime*rate) + usingFraction*value;
		else
			*pConc = value;
	}
	else
	{	// FLIP update
		*pConc += deltime*rate;
	}
#else
	*pConc += deltime*rate;
#endif
	
	return nextTask;
}

#pragma mark ACCESSORS

// Return name of this task
const char *DiffusionTask::TaskName(void)
{	return "diffusion calculations";
}

// Return style of this diffusion this task
const char *DiffusionTask::StyleName(void)
{	switch(style)
	{	case POROELASTICITY_DIFFUSION:
			return "pore pressure";
		case MOISTURE_DIFFUSION:
			return "concentration";
		default:
			break;
	}
	return "unknown diffusion style";
}

// adjust time for given cell size if needed (only for diffusion, other diffusions override)
TransportTask *DiffusionTask::TransportTimeStepFactor(int matid,double *diffCon)
{	*diffCon = theMaterials[matid]->MaximumDiffusion()/(theMaterials[matid]->GetDiffusionCT());
    return nextTask;
}

// return first boundary condition
NodalValueBC *DiffusionTask::GetFirstBCPtr(void) const
{   return firstDiffBC[style];
}
MatPtLoadBC *DiffusionTask::GetFirstFluxBCPtr(void) const
{   return firstDiffFluxBC[style];
}

// return pointer on node to transport field
TransportField *DiffusionTask::GetTransportFieldPtr(NodalPoint *ndpt) const
{   return &(ndpt->gDiff[number]);
}

// particle values and pointers
double DiffusionTask::GetParticleValue(MPMBase *mptr) const
{	return mptr->pDiff[number]->conc;
}
double *DiffusionTask::GetParticleValuePtr(MPMBase *mptr) const
{	return &(mptr->pDiff[number]->conc);
}
double DiffusionTask::GetPrevParticleValue(MPMBase *mptr) const
{	return mptr->pDiff[number]->prevConc;
}
double *DiffusionTask::GetPrevParticleValuePtr(MPMBase *mptr) const
{	return &(mptr->pDiff[number]->prevConc);
}

// Current difference between previous particle value and reference
double DiffusionTask::GetDeltaConcentration(MPMBase *mptr) const { return GetPrevParticleValue(mptr)-reference; }

int DiffusionTask::GetNumber(void) const { return number; }
void DiffusionTask::SetNumber(int taskNum) { number = taskNum; }

int DiffusionTask::GetStyle(void) const { return style; }

#pragma mark STATIC_METHODS

// to check on diffusion or poroelasticity
bool DiffusionTask::HasDiffusion(void)
{	if(diffusion==NULL) return false;
	return diffusion->style==MOISTURE_DIFFUSION;
}
#ifdef POROELASTICITY
bool DiffusionTask::HasPoroelasticity(void)
{	if(diffusion==NULL) return false;
	return diffusion->style==POROELASTICITY_DIFFUSION;
}
#endif
bool DiffusionTask::HasFluidTransport(void) { return diffusion!=NULL; }

// convert poroelasticity MPa to Pa but no change to concentration potentil
double DiffusionTask::RescalePotential(int phaseStyle)
{	if(diffusion==NULL) return 1.;
	return diffusion->style==POROELASTICITY_DIFFUSION ? UnitsController::Scaling(1.e6) : 1 ;
}

// convert Legacy poroelasticity (dV/V)/time to (dV/V)/time (factor=1) and
// convert concetration flux (kg/(m^2-sec) to Legacy (g/(mm^2 sec))
double DiffusionTask::RescaleFlux(void)
{	if(diffusion==NULL) return 1.;
	return diffusion->style==MOISTURE_DIFFUSION ? UnitsController::Scaling(1.e-3) : 1. ;
}

// find diffusion task for a style or NULL is none available
DiffusionTask *DiffusionTask::FindDiffusionTask(int findStyle)
{
	// just first one (and may be NULL)
	if(findStyle<=POROELASTICITY_DIFFUSION)
    {   if(diffusion==NULL) return NULL;
        return diffusion->style==findStyle ? diffusion : NULL;
    }
	
	// look through other diffusions
	DiffusionTask *nextTask = otherDiffusion;
	while(nextTask!=NULL)
	{	if(nextTask->style==findStyle) break;
		nextTask = (DiffusionTask *)nextTask->GetNextTransportTask();
	}
	return nextTask;		// returns NULL if no matches
}

// find number for diffusion style (or -1 if not found)
int DiffusionTask::FindDiffusionNumber(int findStyle)
{
	// just first one (and may be NULL)
	if(findStyle<=POROELASTICITY_DIFFUSION)
	{   if(diffusion==NULL) return -1;
		return diffusion->style==findStyle ? 0 : -1;
	}
	
	// look through other diffusions
	DiffusionTask *nextTask = otherDiffusion;
	while(nextTask!=NULL)
	{	if(nextTask->style==findStyle) break;
		nextTask = (DiffusionTask *)nextTask->GetNextTransportTask();
	}
	return nextTask!=NULL ? nextTask->GetNumber() : -1 ;
}

// Count diffusion tasks. Done when SetConcentration called in MPMBase
// constructor. When counted, set each task's number for later use.
// Subsequent calls do nothing, which implies numDiffusion is set.
void DiffusionTask::CountDiffusionTasks(void)
{	// only need to count once
	if(numDiffusion>=0) return;
	
	// is there one for moisture or poroelasticity?
	numDiffusion = diffusion==NULL ? 0 : 1;
	
	// count other diffusion-like tasks
	DiffusionTask *nextTask = otherDiffusion;
	while(nextTask!=NULL)
	{	nextTask->SetNumber(numDiffusion);
		numDiffusion++;
		nextTask = (DiffusionTask *)nextTask->GetNextTransportTask();
	}
}

// Find the minimum time step of all diffusion tasks
double DiffusionTask::GetMinimumTimeStep(void)
{	// main diffusion task
	double minStep = diffusion!=NULL ? diffusion->GetTimeStep() : 1.e30 ;
	
	// look through other diffusions
	DiffusionTask *nextTask = otherDiffusion;
	while(nextTask!=NULL)
	{	minStep = fmin(minStep,nextTask->GetTimeStep());
		nextTask = (DiffusionTask *)nextTask->GetNextTransportTask();
	}
	return minStep;
}

// Set all diffusion tasks to stop FMPM
// Called only at start of PeriodicXPIC tasks; those needed FMPM are turned on later
void DiffusionTask::SetDiffusionXPIC(bool setting)
{	// main diffusion task
	if(diffusion!=NULL) diffusion->SetUsingTransportXPIC(setting,1.);
	
	// look through other diffusions
	DiffusionTask *nextTask = otherDiffusion;
	while(nextTask!=NULL)
	{	nextTask->SetUsingTransportXPIC(setting,1.);
		nextTask = (DiffusionTask *)nextTask->GetNextTransportTask();
	}
}


