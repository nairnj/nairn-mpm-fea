/******************************************************************************** 
    TransportTask.cpp
    nairn-mpm-fea
    
    Created by John Nairn on 7/18/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
 
 	Transport Calculations with FMPM Option (note FMPM-only version does FLIP too)
   ---------------------------------------------------------------------------------
 	Prepare for Calculations
 		Verify time step for all cells and transport properties (TransportTimeStepFactor())
 		Call Initialize() - print details, and other needs (conduction allocates pTemp on particles
 	Initialization:
 		Set gTValue, gtValueForGrad, gVCT, and gQ on node to zero (in gDiff[] and gCond)
 	Mass and Momentum Task
 		Extrapolate gTValue, gVCT (Task1Extrapolation())
 			If contact implemented, add extrapolations m... and c...
		In reduction, copy gTValue, and gVCT from ghost to real (Task1Reduction())
 			If contact implemented, add extrapolations m... and c...
		(N)ote that materials where csat varies with position extrapolate gTValueRel
			here and implement like gTValue in next task and nodal calculations below - only diffusion)
	Post M&M Extrapolation Task
		Divide gTValue by gVCT on each node (in standard GetTransportNodalValues())
 		Call TransportBCsAndGradients()
 			Copy no BC value and impose grid Tvalue BCs (ImposeValueBCs(time,true))
 			Find grad Tvalue on particles (GetGradients(), ZeroTransportGradients(), AddTransportGradients())
 				If contact implemented, get gradients for m... and c...
		Paste back initial gTValue
 	Update strains first (if used)
 		Read changes in value to input to constitutive laws
 	Grid Forces Task
 		Add transport force to gQ (AddForces())
 			If contact implemented, add for m... and c...
 			Material point classes must return force given gradient and gradient shape functions
 		In reduction, copy gQ from ghost to real (Task1Reduction())
 			If contact implemented, add extrapolations m... and c...
	Post Forces Task
 		Old Method
 			Impose force (grid BC in gQ) and add flux BCs (SetTransportForceAndFluxBCs())
 				If contact implemented, add for m... and c... (only in force part)
 		FMPM Method with FLIP
 			Impose only flux BCs in gQ (SetTransportFluxBCs())
			In transport only, this is done in GridForcesTask.cpp and PostForcesTask.cpp is not used
 	Update Momenta Task
 		Divide gQ by gVCT to get transport rates (UpdateTransport()) and update nodal value
 			If contact implemented, get rates in m... and c...
 		FMPM (with FLIP) approach only
 			Impose grid BCs (TransportGridBCs(mtime,timestep,UPDATE MOMENTUM CALL)
                Calls ImposeValueGridBCs(UPDATE MOMENTUM CALL))
 				For FMPM this is lumped mass matrix preliminary calculations
 	XPIC/FMPM Calculations
 		Calculate gTstar = Sum gTk starting with gTk = m*gTValue (* to format these comments better)
 	Update Particles Task
 		If FMPM with k>1
            Copy gTstar (theta(k)) to gTValue
            Save lumped value in gTstar to implement blending FLIP/FPMP(k>1)
 			Impose grid BCs (TransportGridBCs(mtime,timestep,UPDATE GRID STRAINS CALL)
                Calls ImposeValueGridBCs(UPDATE GRID STRAINS CALL))
 		Each particle: zero rate, value, and pic value, and then extrapolate all from nodes to particle
            value - from theta(k) (if any FMPM(k>1)) or gets pic value if all FLIP (or if FMPM(1))
            If any FPMPM
                If blended FMPM/FLIP, get rate for FLIP and pic (only needed if k>1)
            else
 `              get rate for FLIP
 		Find d(value) from change in extrapolated value (GetDeltaValue())
 		Update particle value using rate (MoveTransportValue())
 		Special Tasks: conduction does adiabatic heating, energy, entropy, and phase transitions
 			poroelasticity does adiabatic pressure
 		Store d(value) on particle for input to next strain updates
 		Find grad Tvalue on particles (GetGradients(), ZeroTransportGradients(), AddTransportGradients())
 			If contact implemented, get for m... and c...
 	Update strains last (if used and if second extrapolation)
 		Read changes in value to input to constitutive laws

	Material Point class support
		Place to store pTValue, pPreviousTValue, and extraolated gradient
		Allocate memory for extrapolated gradient (1 or 2 more gradients if contact implemented)
	Nodal Point class support
		Place to store gTValue, gVCT, and gQ for value, transport "mass", and transport flow rate
		Place to store XPIC extrapolations
		Crack velocity field create structure for contact extrapolations (if needed)
		Material velocity field create structure for contact extrapolations (if needed)
	Material class support
		Include transport properties for k and kappa
		Return "max diffusivity" = max(kappa)/k for use in time step calculations
	Boundary conditions
		Support nodal value and flux conditions
	Contact
		Need material and crack contact methods to adjust transport
********************************************************************************/

#include "stdafx.h"
#include "Custom_Tasks/TransportTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Elements/ElementBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Exceptions/CommonException.hpp"
#include "Boundary_Conditions/NodalValueBC.hpp"
#include "Boundary_Conditions/MatPtLoadBC.hpp"
#include "NairnMPM_Class/XPICExtrapolationTaskTO.hpp"

bool TransportTask::hasXPICOption = false;
int TransportTask::XPICOrder = 0;

// Task lists
TransportTask *transportTasks = NULL;
int numTransport = 0;

// true if either conduction or diffusion are doing contact calculations
bool TransportTask::hasContactEnabled = false;

#pragma mark INITIALIZE

// Constructors
TransportTask::TransportTask()
{	nextTask=NULL;
	transportTimeStep = 1.e30;
	usingXPIC = false;
	usingFraction = 1.;
	doCrelExtrapolation = false;
}

// Destructor (and it is virtual)
TransportTask::~TransportTask() { }

// find time step
void TransportTask::CheckTimeStep(double tst)
{	if(tst<transportTimeStep)
		transportTimeStep = tst;
}

// return time step
double TransportTask::GetTimeStep(void) const { return transportTimeStep; }

#pragma mark MASS AND MOMENTUM EXTRAPOLATIONS

// Task 1 Extrapolation of transport value to the grid for contact (overridden by contact classes)
void TransportTask::Task1ContactExtrapolation(NodalPoint *ndpt,short vfld,int matfld,double argT,double mTpShape) {}

// Task 1 reduction of ghost node to real node for temperature on the grid
TransportTask *TransportTask::Task1Reduction(NodalPoint *real,NodalPoint *ghost)
{	TransportField *gReal = GetTransportFieldPtr(real);
	TransportField *gGhost = GetTransportFieldPtr(ghost);
	gReal->gTValue += gGhost->gTValue;
	gReal->gVCT += gGhost->gVCT;
	if(doCrelExtrapolation) gReal->gTValueRel += gGhost->gTValueRel;
	Task1ContactReduction(real,ghost);
	return nextTask;
}

// Task 1 reduction of transport value to the grid for contact (overridden by contact classes)
void TransportTask::Task1ContactReduction(NodalPoint *real,NodalPoint *ghost) {}

// if active, call transport task to get nodal value
TransportTask *TransportTask::GetTransportNodalValue(NodalPoint *ndptr)
{	if(ndptr->NodeHasNonrigidParticles())
	{	TransportField *gTrans = GetTransportFieldPtr(ndptr);
		gTrans->gTValue /= gTrans->gVCT;
		if(doCrelExtrapolation) gTrans->gTValueRel /= gTrans->gVCT;
		GetContactNodalValue(ndptr);
	}
	return nextTask;
}

// if active, call transport task to get nodal value
void TransportTask::GetContactNodalValue(NodalPoint *ndptr) {}

// Task 1b - impose grid-based transport value BCs
// Changes here may need to added to DiffusionTasks and ConductionTask(s) as well
void TransportTask::ImposeValueBCs(double stepTime,bool copyFirst)
{   // exit if no BCs
 	NodalValueBC *nextBC = GetFirstBCPtr();
    if(nextBC==NULL) return;
    
    // Copy no-BC transport value
    int i;
	if(copyFirst)
	{	while(nextBC!=NULL)
		{   i = nextBC->GetNodeNum(stepTime);
			if(i!=0) nextBC->CopyNodalValue(nd[i],GetTransportFieldPtr(nd[i]));
			nextBC = (NodalValueBC *)nextBC->GetNextObject();
		}
        nextBC = GetFirstBCPtr();
	}
    
    // Zero them all
    while(nextBC!=NULL)
    {   i = nextBC->GetNodeNum(stepTime);
        if(i!=0)
        {   TransportField *gTrans = GetTransportFieldPtr(nd[i]);
            gTrans->gTValue = 0.;
        }
        nextBC = (NodalValueBC *)nextBC->GetNextObject();
    }
    
    // Now add all transport values to nodes with value BCs
    nextBC = GetFirstBCPtr();
    while(nextBC!=NULL)
    {   i = nextBC->GetNodeNum(stepTime);
        if(i!=0)
        {   TransportField *gTrans = GetTransportFieldPtr(nd[i]);
            gTrans->gTValue += nextBC->BCValue(stepTime);
        }
        nextBC = (NodalValueBC *)nextBC->GetNextObject();
    }
}

// Paste back value that were recently copied
TransportTask *TransportTask::RestoreValueBCs(void)
{   // exit if no BCs
    NodalValueBC *nextBC = GetFirstBCPtr();
    if(nextBC==NULL) return nextTask;
	
	// Paste back noBC transport value
    int i;
	while(nextBC!=NULL)
	{	i = nextBC->GetNodeNum(mtime);
		if(i!=0) nextBC->PasteNodalValue(nd[i],GetTransportFieldPtr(nd[i]));
		nextBC = (NodalValueBC *)nextBC->GetNextObject();
	}
	
	return nextTask;
}

// Task 1b - get gradients in transport value on particles
// throws CommonException()
TransportTask *TransportTask::GetGradients(double stepTime)
{
	CommonException *transErr = NULL;
#ifdef CONST_ARRAYS
	int ndsArray[MAX_SHAPE_NODES];
	double fn[MAX_SHAPE_NODES],xDeriv[MAX_SHAPE_NODES],yDeriv[MAX_SHAPE_NODES],zDeriv[MAX_SHAPE_NODES];
#else
	int ndsArray[maxShapeNodes];
	double fn[maxShapeNodes],xDeriv[maxShapeNodes],yDeriv[maxShapeNodes],zDeriv[maxShapeNodes];
#endif
	
	// in case 2D planar
	for(int i=0;i<maxShapeNodes;i++) zDeriv[i] = 0.;
	
	// Find gradients on the nonrigid particles
#pragma omp parallel for private(ndsArray,fn,xDeriv,yDeriv) firstprivate(zDeriv)
	for(int p=0;p<nmpmsNR;p++)
	{	if(mpm[p]->InReservoir()) continue;
		try
		{   // find shape functions and derviatives
			MPMBase *mptr = mpm[p];
			const ElementBase *elref = theElements[mptr->ElemID()];
			int *nds = ndsArray;
			elref->GetShapeGradients(fn,&nds,xDeriv,yDeriv,zDeriv,mptr);
			int i,numnds = nds[0];
			
			// Zero transport gradients on particles
			ZeroTransportGradients(mptr);
			
			// Find gradients from current temperatures
			Vector deriv;
			for(i=1;i<=numnds;i++)
			{	deriv = MakeVector(xDeriv[i],yDeriv[i],zDeriv[i]);
				AddTransportGradients(mptr,&deriv,nd[nds[i]],mptr->vfld[i]);
			}
		}
		catch(CommonException& err)
		{   if(transErr!=NULL)
			{
#pragma omp critical (error)
				transErr = new CommonException(err);
			}
		}
	}
	
	// throw any errors
	if(transErr!=NULL) throw *transErr;
	return nextTask;
}

// Zero gradients on the particles for contact calculations (overridden by contact classes)
void TransportTask::ZeroTransportContactGradients(MPMBase *mptr) {}

// Add gradients on the particles (overridden by contact classes)
void TransportTask::AddTransportContactGradients(MPMBase *mptr,Vector *deriv,NodalPoint *ndptr,short vfld) {}

#pragma mark GRID FORCES EXTRAPOLATIONS

// Reduction of transport forces to the grid
TransportTask *TransportTask::ForcesReduction(NodalPoint *real,NodalPoint *ghost)
{	TransportField *gReal = GetTransportFieldPtr(real);
	TransportField *gGhost = GetTransportFieldPtr(ghost);
	gReal->gQ += gGhost->gQ;
	ForcesContactReduction(real,ghost);
	return nextTask;
}

// Reduction of transport forces to the grid for contact (overridden by contact classes)
void TransportTask::ForcesContactReduction(NodalPoint *real,NodalPoint *ghost) {}

// find forces for conduction calculation (N-mm/sec = mJ/sec) (non-rigid particles only)
void TransportTask::AddContactForces(NodalPoint *ndptr,MPMBase *mptr,double sh,double dshdx,
										 double dshdy,double dshdz,TransportProperties *t,short vfld,int matfld)
{}


// add flux to transport flow force to active nodes only
// postRateCalc means gCond.gQ has been divided by its mass
void TransportTask::AddFluxCondition(NodalPoint *ndptr,double extraFlux,bool postRateCalc)
{	if(ndptr->NodeHasNonrigidParticles())
    {   TransportField *gTrans = GetTransportFieldPtr(ndptr);
		if(postRateCalc) extraFlux /= gTrans->gVCT;
        gTrans->gQ += extraFlux;
    }
}

// Impose grid-based transport value BCs after updating values on the grid in UpdateMomentumTask.
// FLIP: stype==UPDATE_GRID_STRAINS_CALL never called
// FMPM: callStyle==UPDATE_GRID_STRAINS_CALL only called with order>1, no need to update rate
//			QReaction is recalculated
// FLIP/FMPM mix: UPDATE_GRID_STRAINS_CALL only called with order>1, currently skips rate
//			and recalculates QReaction. Not certain yet if this is correct or not
TransportTask *TransportTask::ImposeValueGridBCs(double bctime,double deltime,int callStyle)
{   // Skip if no BCs
    // initialize for transport BC
    NodalValueBC *nextBC = GetFirstBCPtr();
    if(nextBC==NULL) return nextTask;
    
	// --------- set value and consistent rates for grid transport BCs ------------
	int i;
	TransportField *gTrans;
 	
	// The rest is update momentum call - set BC in lumped capacity value on the grid
	// Set rates too
	while(nextBC!=NULL)
	{	i = nextBC->GetNodeNum(mtime);
		if(i!=0)
		{	gTrans = GetTransportFieldPtr(nd[i]);
			gTrans->gFirstBC = true;
            if(callStyle==UPDATE_GRID_STRAINS_CALL)
                nextBC->InitQReaction(1.-usingFraction);
            else
                nextBC->InitQReaction();
		}
		nextBC = (NodalValueBC *)nextBC->GetNextObject();
	}
	
	// Set Value to Sum(bc) Value(BC)
	// Change gQ to (Sum(bc) Value(BC) - Value(0))/dt
	//    But: Value was just changed in momentum update to Value(up) = Value(0) + gQ*dT
	//    Solution: Add (Sum(bc) Value(BC) - Value(up))/dt to gQ to get (Sum(bc) Value(BC) - Value(0))/dt
	//          (Sum(bc) Value(BC) - Value(0) - gQ*dt)/dt + dt =
	// Add gVCT*(Sum(bc) Value(BC) - Value(0))/dt to sum of BC reactions
	double dv;
	nextBC = GetFirstBCPtr();
	while(nextBC!=NULL)
	{	i = nextBC->GetNodeNum(mtime);
		if(i!=0)
		{	gTrans = GetTransportFieldPtr(nd[i]);
			
			// first one on this node
			if(gTrans->gFirstBC)
			{	// temperature rate -T(0)/dt first one only
				dv = -gTrans->gTValue/deltime;
				// Note summing is correct because done after momentum update (see above)
                if(callStyle==UPDATE_GRID_STRAINS_CALL)
                    nextBC->SuperposeQReaction(usingFraction*gTrans->gVCT*dv);
                else
                {   gTrans->gQ += dv;
                    nextBC->SuperposeQReaction(gTrans->gVCT*dv);
                }
				
				// clear value for first one
				gTrans->gTValue = 0.;
				
				// only do this once
				gTrans->gFirstBC = false;
			}
			
			// add to concentration, rate, and Qreaction
			double bcT = nextBC->BCValue(bctime);
			gTrans->gTValue += bcT;
			dv = bcT/deltime;
            if(callStyle==UPDATE_GRID_STRAINS_CALL)
                nextBC->SuperposeQReaction(usingFraction*gTrans->gVCT*dv);
			else
            {   gTrans->gQ += dv;
                nextBC->SuperposeQReaction(gTrans->gVCT*dv);
            }
		}
		nextBC = (NodalValueBC *)nextBC->GetNextObject();
	}
	
	return nextTask;
}

// Add Flux BCs
TransportTask *TransportTask::SetTransportFluxBCs(void)
{
	MatPtLoadBC *nextFlux = GetFirstFluxBCPtr();
	while(nextFlux!=NULL)
		nextFlux = nextFlux->AddMPFluxBC(mtime);
	
	// next task
	return nextTask;
}

#pragma mark UPDATE MOMENTA TASK AND CONTACT FLOW

// Get transport rate but active nodes only
TransportTask *TransportTask::UpdateTransport(NodalPoint *ndptr,double deltime)
{	if(ndptr->NodeHasNonrigidParticles())
	{	TransportField *gTrans = GetTransportFieldPtr(ndptr);
		gTrans->gQ /= gTrans->gVCT;
		gTrans->gTValue += gTrans->gQ*timestep;
		TransportContactRates(ndptr,deltime);
	}
	return nextTask;
}

// Reduction of transport forces to the grid for contact (overridden by contact classes)
void TransportTask::TransportContactRates(NodalPoint *ndptr,double deltime) {}

#pragma mark UPDATE PARTICLES TASK

// increment temperature rate on the particle (only used in transport XPIC)
// only implemented for global temperature - need update to handle contact heat flow
TransportTask *TransportTask::InitializeForXPIC(NodalPoint *ndptr,double timestep,int xpicOption) const
{	TransportField *gTrans = GetTransportFieldPtr(ndptr);
	gTrans->gTnext = 0.;
	// skip if ghost node (ghost only needs gTnext=0)
	if(timestep>0.)
	{	double m = TransportTask::XPICOrder;
		gTrans->gTprev = m*gTrans->gTValue;
		gTrans->gTstar = gTrans->gTprev;
	}
	return nextTask;
}

// increment temperature rate on the particle
double TransportTask::IncrementTransportRate(NodalPoint *ndptr,double shape,short vfld,int matfld) const
{	TransportField *gTrans = GetTransportFieldPtr(ndptr);
	return gTrans->gQ*shape;
}

// after extrapolated, give task to change rate and value if needed
void TransportTask::AdjustRateAndValue(MPMBase *mptr,double &value,
									   double &rate,double &lumpedValue,double deltime) const {}

// increment particle concentration (time is always timestep)
TransportTask *TransportTask::MoveTransportValue(MPMBase *mptr,double deltime,double rate,double value) const
{   double *pValue = GetParticleValuePtr(mptr);
	if(usingXPIC)
	{	// FMPM update just replace the value or blend with FLIP
		if(usingFraction<1.)
			*pValue = (1.-usingFraction)*(*pValue+deltime*rate) + usingFraction*value;
		else
			*pValue = value;
	}
	else
	{	// FLIP update
		*pValue += deltime*rate;
	}
    return nextTask;
}

#pragma mark UPDATE PARTICLE STRAIN TASK

// return increment transport value
double TransportTask::IncrementValueExtrap(NodalPoint *ndptr,double shape,short vfld,int matfld) const
{	TransportField *gTrans = GetTransportFieldPtr(ndptr);
	return gTrans->gTValue*shape;
}

// return increment transport lumped value (only call when using FMPM with k>1
double TransportTask::IncrementLumpedValueExtrap(NodalPoint *ndptr,double shape,short vfld,int matfld) const
{	TransportField *gTrans = GetTransportFieldPtr(ndptr);
	return gTrans->gTstar*shape;
}

// after extrapolated, find change this update on particle extrapolated from the grid
// Note that all "otherDiffusion" tasks must override to not set *dV (it is only
//		used for dT, dc, and dp) and to store delta value in history variable
//		if needed. This need is when it does not just return dV value.
void TransportTask::GetDeltaValue(MPMBase *mptr,double pValueExtrap,double *dV) const
{   double *pPrevValue = GetPrevParticleValuePtr(mptr);
    double dValue = pValueExtrap-(*pPrevValue);
    *pPrevValue = pValueExtrap;
    *dV = dValue;
}

#pragma mark ACCESSORS

// Return name of this task
const char *TransportTask::TaskName(void) { return "transport calculations"; }

// get the next task
TransportTask *TransportTask::GetNextTransportTask(void) const { return nextTask; }

// only used by diffusion tasks (here for access in TransportTask methods)
int TransportTask::GetNumber(void) const { return -1; }
void TransportTask::SetNumber(int taskNum) { }

// to activate XPIC
void TransportTask::SetUsingTransportXPIC(bool setting,double fractionFMPM)
{	usingXPIC = setting;
	usingFraction = fractionFMPM;
}
bool TransportTask::IsUsingTransportXPIC(void) const { return usingXPIC; }
bool TransportTask::IsUsingTransportXPIC(double &fraction) const
{	fraction = usingFraction;
	return usingXPIC;
}
// if blended FLIP/FMPM(k>1) and want grid extrapolation blended, then
//   return true and set fraction, otherwise return false
// Other alternatives
//	always return false: uses FMPM(k) extrapolation when available
//  set fraction=0 when returning true: uses lumped values
bool TransportTask::ShouldBlendFromGrid(double &fraction) const
{
	// if FLIP or pure FMPM(k>0) or blended FLIP/FMPM(1) then
	//   nothing available to blend so return false
	if(!usingXPIC) return false;
	if(usingFraction>=1.) return false;
	if(XPICTransportTask->GetXPICOrder()<2) return false;
	
	// Here for blended FLIP/FMPM(k>1), set fraction and return
	fraction = usingFraction;
	return true;
}

#pragma mark CLASS METHODS

// get transport values on nodes in post mass and momentum extrapolation tasks
void TransportTask::GetTransportValues(NodalPoint *ndptr)
{
	TransportTask *nextTransport=transportTasks;
	while (nextTransport != NULL)
		nextTransport = nextTransport->GetTransportNodalValue(ndptr);
}

// Impose transport BCs and extrapolate gradients to the particles
void TransportTask::TransportBCsAndGradients(double bctime)
{
	// impose BCs in lumped capacity values
	// and get gradients from lumped capacity values
	TransportTask *nextTransport=transportTasks;
	while(nextTransport!=NULL)
	{   nextTransport->ImposeValueBCs(bctime,true);
		nextTransport = nextTransport->GetGradients(bctime);
	}

	// restore the values on BC nodes that were changed above
	nextTransport=transportTasks;
	while(nextTransport!=NULL)
		nextTransport = nextTransport->RestoreValueBCs();
}

// Set transport BCs (not parallel because small and possible use of function/global variables)
void TransportTask::TransportForceBCs(double dtime)
{   TransportTask *nextTransport=transportTasks;
	while(nextTransport!=NULL)
		nextTransport=nextTransport->SetTransportFluxBCs();
}

// Called suring the momentum update
// Get grid transport rates and updated value (lumped mass matrix method)
// For contact, do material and crack contact updates too
void TransportTask::UpdateTransportOnGrid(NodalPoint *ndptr)
{
	TransportTask *nextTransport=transportTasks;
	while(nextTransport!=NULL)
		nextTransport = nextTransport->UpdateTransport(ndptr,timestep);
}

// Impose transport BCs on the grid. It is called in UpdateMomentaTask
//     to use the new grid values
// If FMPM activated and order>1, it is called again with v* on grid from
//    XPIC tasks. This task imposes BCs in v*. No change is needed for
//    transport tasks not using FMPM (because they already had BCs imposed
//	  in lumped grid velocity)
// Finding rate might be an issue for blended FLIP/FMPM?
void TransportTask::TransportGridBCs(double bctime,double deltime,int callStyle)
{
	TransportTask *nextTransport;
	
	if(callStyle==UPDATE_MOMENTUM_CALL)
	{	nextTransport=transportTasks;
		while(nextTransport!=NULL)
		{	if(callStyle==UPDATE_GRID_STRAINS_CALL)
			{	// only called when FMPM activated, but only call if using FMPM
				if(nextTransport->IsUsingTransportXPIC())
				{	// imposes BCs
					nextTransport->ImposeValueGridBCs(bctime,deltime,callStyle);
				}
			}
			else
			{	// always call in momentum task
				nextTransport->ImposeValueGridBCs(bctime,deltime,callStyle);
			}
			nextTransport = nextTransport->GetNextTransportTask();
		}
	}
}
