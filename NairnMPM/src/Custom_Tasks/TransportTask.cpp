/******************************************************************************** 
    TransportTask.cpp
    nairn-mpm-fea
    
    Created by John Nairn on 7/18/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
 
	Transport calculations
   -------------------------
    Prepare for Calculations
		Verify time step for all cells and transport properties (TransportTimeStepFactor())
		Call Initialize() - allocate particle gradients, print details, and other needs
	Initialization:
		Set gTValue, gMTp, and gQ on node to zero;
	Mass and Momentum Task
		Extrapolate gTValue, gMTp (Task1Extrapolation())
			If contact implemented, add extrapolations m... and c...
		In reduction, copy gTValue and gMTp from ghost to real (Task1Reduction())
			If contact implemented, add extrapolations m... and c...
		Post Extrapolation Task
			Divide gTvale by gMTp (GetTransportNodalValue())
				If contact implemented, get values m... and c...
			Impose grid Tvalue BCs (ImposeValueBCs())
				If contact implemented, apply to m... and c...
			Find grad Tvalue on particles (GetGradients(), ZeroTransportGradients(), AddTransportGradients())
				If contact implemented, get for m... and c...
	Update strains first (if used)
		See update strains on particles below
 	Grid Forces Task
		Add transport force to gQ (AddForces())
			If contact implemented, add for m... and c...
		Material point classes must return force give gradient and gradient shape functinos
		In reduction, copy gQ from ghost to real (Task1Reduction())
			If contact implemented, add extrapolations m... and c...
		Post Forces Task
			Finish grid BCs and impose flux BCs in gQ (SetTransportForceBCs())
				If contact implemented, impose in m... and c...
 	Update Momenta Task
		Divide gQ by gMtp to get transport rates (GetTransportRates())
			If contact implemented, get rates in m... and c...
	Update Particles Task
        Each particle: zero rate, extrapolate from nodes to particle, then
            update particle (IncrementTransportRate(),MoveTransportRate()).
        Add adiabatic term to particle transport value
	Update strains last (if used and if second extrapolation)
		Update transport value (UpdateNodalValues())
            If contact implemented, get values m... and c...
		See update strains on particles below
 
	Update strains on particles
		In UpdateStrain() for material point classes, extrapolate grid, cvf, or mvf
			to particle using IncrementValueExtrap()
		Find dT = extrapolated value minus previous extrapolated value and then
			store bew extrapolated value in previous value (GetDeltaValue())
		This task is done because contititutive laws work better when transport value comes
			from grid extrapolation instead of particle value. The current grid
			based result is stored in the previous value. The particle
			value, however, is the one used in transport calculations. Using the
			grid based one causes numerical diffision.
 
	Material Point class support
		Place to store pTValue, pPreviousTValue, and extraolated gradient
		Allocate memory for extrapolated gradient (1 or 2 more gradients if contact implemented)
	Nodal Point class support
		Place to store gTValue, gMTp, and gQ for value, transport "mass", and transport flow rate
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

TransportTask *transportTasks=NULL;

// true if either conduction or diffusion are doing contact calculations
bool TransportTask::hasContactEnabled=false;

#pragma mark INITIALIZE

// Constructors
TransportTask::TransportTask()
{	nextTask=NULL;
}

// Destructor (and it is virtual)
TransportTask::~TransportTask() { }

#pragma mark MASS AND MOMENTUM EXTRAPOLATIONS

// Task 1 Extrapolation of temperature to the grid
// Only called for non-rigid materials
TransportTask *TransportTask::Task1Extrapolation(NodalPoint *ndpt,MPMBase *mptr,double shape,short vfld,int matfld)
{	double pTValue;
	double mTpShape = GetTransportMassAndValue(mptr,&pTValue)*shape;
	double mTpTValueShape = pTValue*mTpShape;
	TransportField *gTrans = GetTransportFieldPtr(ndpt);
	gTrans->gTValue += mTpTValueShape;
	gTrans->gMTp += mTpShape;
	Task1ContactExtrapolation(ndpt,vfld,matfld,mTpTValueShape,mTpShape);
	return nextTask;
}

// Task 1 Extrapolation of transport value to the grid for contact (overridden by contact classes)
void TransportTask::Task1ContactExtrapolation(NodalPoint *ndpt,short vfld,int matfld,double argT,double mTpShape) {}

// Task 1 reduction of ghost node to real node for temperature on the grid
TransportTask *TransportTask::Task1Reduction(NodalPoint *real,NodalPoint *ghost)
{	TransportField *gReal = GetTransportFieldPtr(real);
	TransportField *gGhost = GetTransportFieldPtr(ghost);
	gReal->gTValue += gGhost->gTValue;
	gReal->gMTp += gGhost->gMTp;
	Task1ContactReduction(real,ghost);
	return nextTask;
}

// Task 1 reduction of transport value to the grid for contact (overridden by contact classes)
void TransportTask::Task1ContactReduction(NodalPoint *real,NodalPoint *ghost) {}

// if active, call transport task to get nodal value
TransportTask *TransportTask::GetTransportNodalValue(NodalPoint *ndptr)
{	if(ndptr->NodeHasNonrigidParticles())
	{	TransportField *gTrans = GetTransportFieldPtr(ndptr);
		gTrans->gTValue /= gTrans->gMTp;
		GetContactNodalValue(ndptr);
	}
	return nextTask;
}

// if active, call transport task to get nodal value
void TransportTask::GetContactNodalValue(NodalPoint *ndptr) {}

// Task 1b - impose grid-based transport value BCs
void TransportTask::ImposeValueBCs(double stepTime)
{
    int i;
    
    // Copy no-BC transport value
    NodalValueBC *nextBC = GetFirstBCPtr();
    while(nextBC!=NULL)
    {   i = nextBC->GetNodeNum(stepTime);
        if(i!=0) nextBC->CopyNodalValue(nd[i]);
        nextBC = (NodalValueBC *)nextBC->GetNextObject();
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
        if(postRateCalc) extraFlux /= gTrans->gMTp;
        gTrans->gQ += extraFlux;
    }
}

// adjust forces at grid points with concentration BCs to have rates be correct
// to carry extrapolated concentrations (before impose BCs) to the correct
// one selected by grid based BC
TransportTask *TransportTask::SetTransportForceBCs(double deltime)
{
	int i;
	TransportField *gTrans;
	
	// --------- consistent forces for grid concentration BCs ------------
	
	// Paste back noBC concentration
	NodalValueBC *nextBC = GetFirstBCPtr();
	while(nextBC!=NULL)
	{	i = nextBC->GetNodeNum(mtime);
		if(i!=0) nextBC->PasteNodalValue(nd[i]);
		nextBC = (NodalValueBC *)nextBC->GetNextObject();
	}
	
	// Set force to - VC(no BC)/timestep
	nextBC = GetFirstBCPtr();
	while(nextBC!=NULL)
	{	i = nextBC->GetNodeNum(mtime);
		if(i!=0)
		{	gTrans = GetTransportFieldPtr(nd[i]);
			gTrans->gQ = -gTrans->gMTp*gTrans->gTValue/deltime;
		}
		nextBC = (NodalValueBC *)nextBC->GetNextObject();
	}
	
	// Now add each superposed concentration (* volume) BC at incremented time
	nextBC = GetFirstBCPtr();
	while(nextBC!=NULL)
	{	i = nextBC->GetNodeNum(mtime);
		if(i!=0)
		{	gTrans = GetTransportFieldPtr(nd[i]);
			gTrans->gQ += gTrans->gMTp*nextBC->BCValue(mtime)/deltime;
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

#ifdef CONTACT_HEAT_FLOW
// adjust for material contact
TransportTask *TransportTask::MatContactFlowCalculations(MatVelocityField *mvf,NodalPoint *ndptr,
														 CrackVelocityField *cvf,double surfaceArea,bool inContact)
{	return nextTask;
}

// adjust for material contact
TransportTask *TransportTask::CrackContactFlowCalculations(CrackVelocityField *mvf,NodalPoint *ndptr,double surfaceArea,int inContact)
{	return nextTask;
}
#endif

// Get transport rate but active nodes only
TransportTask *TransportTask::GetTransportRates(NodalPoint *ndptr,double deltime)
{	if(ndptr->NodeHasNonrigidParticles())
	{	TransportField *gTrans = GetTransportFieldPtr(ndptr);
		gTrans->gQ /= gTrans->gMTp;
		TransportContactRates(ndptr,deltime);
	}
	return nextTask;
}

// Reduction of transport forces to the grid for contact (overridden by contact classes)
void TransportTask::TransportContactRates(NodalPoint *ndptr,double deltime) {}

#pragma mark UPDATE PARTICLES TASK

// increment temperature rate on the particle
TransportTask *TransportTask::IncrementTransportRate(NodalPoint *ndptr,double shape,double &rate,short vfld,int matfld) const
{	TransportField *gTrans = GetTransportFieldPtr(ndptr);
	rate += gTrans->gQ*shape;
	return nextTask;
}

// increment particle concentration (time is always timestep)
TransportTask *TransportTask::MoveTransportValue(MPMBase *mptr,double deltime,double rate) const
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

#pragma mark UPDATE STRAIN LAST TASK

// if needed for SZS or USAVG, update temperature on the grid (tempTime is always timestep)
TransportTask *TransportTask::UpdateNodalValues(double tempTime)
{	// add for each node
	for(int i=1;i<=nnodes;i++)
	{   if(nd[i]->NodeHasNonrigidParticles())
		{	TransportField *gTrans = GetTransportFieldPtr(nd[i]);
			gTrans->gTValue += gTrans->gQ*tempTime;
			UpdateContactNodalValues(nd[i],tempTime);
		}
	}
	return nextTask;
}

// Reduction of transport forces to the grid for contact (overridden by contact classes)
void TransportTask::UpdateContactNodalValues(NodalPoint *ndptr,double tempTime) {}

#pragma mark ACCESSORS

// Return name of this task
const char *TransportTask::TaskName(void) { return "transport calculations"; }

// get the next task
TransportTask *TransportTask::GetNextTransportTask(void) const { return nextTask; }

