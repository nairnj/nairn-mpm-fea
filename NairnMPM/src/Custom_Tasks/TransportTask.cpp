/******************************************************************************** 
    TransportTask.cpp
    nairn-mpm-fea
    
    Created by John Nairn on 7/18/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
 
 	Transport Calculations with FMPM Option
   ------------------------------------------
 	Prepare for Calculations
 		Verify time step for all cells and transport properties (TransportTimeStepFactor())
 		Call Initialize() - allocate tranport data (pTemp and pDiffusion), print details, and other needs
 	Initialization:
 		Set gTValue, gVCT, and gQ on node to zero (in gDiff and gCond)
 	Mass and Momentum Task
 		Extrapolate gTValue, gVCT (Task1Extrapolation())
 			If contact implemented, add extrapolations m... and c...
		In reduction, copy gTValue and gVCT from ghost to real (Task1Reduction())
 			If contact implemented, add extrapolations m... and c...
	Post M&M Extrapolation Task
		Divide gTValue by gVCT (GetTransportValues() -> GetTransportNodalValues())
			If contact implemented, get values m... and c...
 		Call TransportBCsAndGradients()
 			Copy no BC value and impose grid Tvalue BCs (ImposeValueBCs())
 				If contact implemented, apply to m... and c...
 			Find grad Tvalue on particles (GetGradients(), ZeroTransportGradients(), AddTransportGradients())
 				If contact implemented, get gradients for m... and c...
		Paste back initial gTValue (but only done in FMPM version, even if using FLIP it that version?)
 	Update strains first (if used)
 		Read changes in value to input to constitutive laws
 	Grid Forces Task
 		Add transport force to gQ (AddForces())
 			If contact implemented, add for m... and c...
 			Material point classes must return force given gradient and gradient shape functinos
 		In reduction, copy gQ from ghost to real (Task1Reduction())
 			If contact implemented, add extrapolations m... and c...
	Post Forces Task
 		Old Method
 			Impose force (grid BC in gQ) and add flux BCs (SetTransportForceAndFluxBCs())
 				If contact implemented, add for m... and c... (only in force part)
 		FMPM Method
 			Impose only flux BCs in gQ (SetTransportFluxBCs())
 	Update Momenta Task
 		Divide gQ by gVCT to get transport rates (UpdateTransport()) and update nodal value
 			If contact implemented, get rates in m... and c...
 		FMPM approach onlu
 			Impose grid BCs (ImposeValueGridBCs(UPDATE_MOMENTUM_CALL))
 				This is lumped mass matrix preliminary calculations
 	XPIC/FMPM Calculations
 		Calculate gTstar = Sum gTk starting with gTk = m*gTValue
 	Update Particles Task
 		Copy gTstar to gTValue (if order>1)
 			Impose grid BCs (ImposeValueBCs(time,false)) (optional)
 		Each particle: zero rate and value, extrapolate both from nodes to particle
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

// Task list
TransportTask *transportTasks=NULL;

// true if either conduction or diffusion are doing contact calculations
bool TransportTask::hasContactEnabled = false;

#pragma mark INITIALIZE

// Constructors
TransportTask::TransportTask()
{	nextTask=NULL;
	transportTimeStep = 1.e30;
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
		GetContactNodalValue(ndptr);
	}
	return nextTask;
}

// if active, call transport task to get nodal value
void TransportTask::GetContactNodalValue(NodalPoint *ndptr) {}

// Task 1b - impose grid-based transport value BCs
void TransportTask::ImposeValueBCs(double stepTime,bool copyFirst)
{
    int i;
	NodalValueBC *nextBC;
    
    // Copy no-BC transport value
	if(copyFirst)
	{	nextBC = GetFirstBCPtr();
		while(nextBC!=NULL)
		{   i = nextBC->GetNodeNum(stepTime);
			if(i!=0) nextBC->CopyNodalValue(nd[i]);
			nextBC = (NodalValueBC *)nextBC->GetNextObject();
		}
	}
    
    // Zero them all
    nextBC = GetFirstBCPtr();
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
	{	try
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

// adjust forces at grid points with concentration BCs to have rates be correct
// to carry extrapolated concentrations (before impose BCs) to the correct
// one selected by grid based BC
TransportTask *TransportTask::SetTransportForceAndFluxBCs(double deltime)
{
	// --------- consistent forces for grid transport BCs ------------
	// note that conducution overrides and gets reaction energy
	int i;
	TransportField *gTrans;

	// Paste back noBC transport value
	NodalValueBC *nextBC = GetFirstBCPtr();
	while(nextBC!=NULL)
	{	i = nextBC->GetNodeNum(mtime);
		if(i!=0) nextBC->PasteNodalValue(nd[i]);
		nextBC = (NodalValueBC *)nextBC->GetNextObject();
	}
	
	// Set force to - Ci*Ti(no BC)/timestep
	nextBC = GetFirstBCPtr();
	while(nextBC!=NULL)
	{	i = nextBC->GetNodeNum(mtime);
		if(i!=0)
		{	gTrans = GetTransportFieldPtr(nd[i]);
			gTrans->gQ = -gTrans->gVCT*gTrans->gTValue/deltime;
		}
		nextBC = (NodalValueBC *)nextBC->GetNextObject();
	}
	
	// Now add each superposed BC (ci*TiBC/timestep) at incremented time
	nextBC = GetFirstBCPtr();
	while(nextBC!=NULL)
	{	i = nextBC->GetNodeNum(mtime);
		if(i!=0)
		{	gTrans = GetTransportFieldPtr(nd[i]);
			gTrans->gQ += gTrans->gVCT*nextBC->BCValue(mtime)/deltime;
		}
		nextBC = (NodalValueBC *)nextBC->GetNextObject();
	}
	
	// --------- concentration flux BCs -------------
	MatPtLoadBC *nextFlux = GetFirstFluxBCPtr();
	while(nextFlux!=NULL)
		nextFlux = nextFlux->AddMPFluxBC(mtime);
	
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

// increment temperature rate on the particle
double TransportTask::IncrementTransportRate(NodalPoint *ndptr,double shape,short vfld,int matfld) const
{	TransportField *gTrans = GetTransportFieldPtr(ndptr);
	return gTrans->gQ*shape;
}

// increment particle concentration (time is always timestep)
TransportTask *TransportTask::MoveTransportValue(MPMBase *mptr,double deltime,double rate,double value) const
{   double *pValue = GetParticleValuePtr(mptr);
	*pValue += deltime*rate;
    return nextTask;
}

#pragma mark UPDATE PARTICLE STRAIN TASK

// return increment transport rate
double TransportTask::IncrementValueExtrap(NodalPoint *ndptr,double shape,short vfld,int matfld) const
{	TransportField *gTrans = GetTransportFieldPtr(ndptr);
	return gTrans->gTValue*shape;
}

// after extrapolated, find change this update on particle
double TransportTask::GetDeltaValue(MPMBase *mptr,double pValueExtrap) const
{   double *pPrevValue = GetPrevParticleValuePtr(mptr);
    double dValue = pValueExtrap-(*pPrevValue);
    *pPrevValue = pValueExtrap;
    return dValue;
}

#pragma mark ACCESSORS

// Return name of this task
const char *TransportTask::TaskName(void) { return "transport calculations"; }

// get the next task
TransportTask *TransportTask::GetNextTransportTask(void) const { return nextTask; }

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
}

// Set transport BCs (not parallel because small and possible use of function/global variables)
void TransportTask::TransportForceBCs(double dtime)
{
    TransportTask *nextTransport=transportTasks;
    while(nextTransport!=NULL)
        nextTransport=nextTransport->SetTransportForceAndFluxBCs(dtime);
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


