/********************************************************************************
	UpdateParticlesTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	The tasks are:
	-------------
	* Extrpolate acceleration, velocity, and spin terms from grid
	  to the particle
	* When activated, extrapolate transport rates to particle
	* Update particle position, velocity, and spin
	* Update particle temperature and concentration
	  (if an adiabatic heating term, add to particle temperature)
	* After main loop, update position of all rigid particles
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/UpdateParticlesTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Custom_Tasks/TransportTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "NairnMPM_Class/XPICExtrapolationTask.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "NairnMPM_Class/XPICExtrapolationTaskTO.hpp"

#pragma mark CONSTRUCTORS

UpdateParticlesTask::UpdateParticlesTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Update particle position, velocity, temp, and conc
// throws CommonException()
bool UpdateParticlesTask::Execute(int taskOption)
{
	CommonException *upErr = NULL;
#ifdef CONST_ARRAYS
	int ndsArray[MAX_SHAPE_NODES];
	double fn[MAX_SHAPE_NODES];
#else
	int ndsArray[maxShapeNodes];
	double fn[maxShapeNodes];
#endif

	// data structure for extrapolations
	GridToParticleExtrap gp;
	
	// Damping terms on the grid or on the particles
	//      particleAlpha   =  pdamping(t)
	//      gridAlpha       = damping(t)
	double particleAlpha = bodyFrc.GetParticleDamping(mtime);
	double gridAlpha = bodyFrc.GetGridDamping(mtime);
	
	// Calculate velocities on the grid
	int m = bodyFrc.GetXPICOrder();
	if(m>1)
	{	// FLIP(k>1) or XPIC(k>1) always has XPICMechanicsTask
		XPICMechanicsTask->Execute(0);
		// velocity BCs (in case used)
		if(bodyFrc.GridBCOption()!=GRIDBC_LUMPED_ONLY)
			NodalVelBC::GridVelocityConditions(UPDATE_GRID_STRAINS_CALL);
	}
	else
	{	// FLIP and PIC (FMPM(1) or XPIC(1)) use lumped mass matrix here
#pragma omp parallel for
		for(int i=1;i<=*nda;i++)
			nd[nda[i]]->GridValueCalculation(VELOCITY_FOR_STRAIN_UPDATE);
	}

	// change sign of m if not FMPM
	if(!bodyFrc.UsingFMPM()) m = -m;
	
#ifdef TRANSPORT_FMPM
	// get grid transport values
	if(XPICTransportTask!=NULL)
	{	if(XPICTransportTask->GetXPICOrder()>1)
		{
#pragma omp parallel for
			for(int i=1;i<=*nda;i++)
				XPICTransportTask->CopyXStar(nd[nda[i]]);
			TransportTask::TransportGridBCs(mtime,timestep,UPDATE_GRID_STRAINS_CALL);
		}
	}
#endif

	// Update particle position, velocity, temp, and conc
#pragma omp parallel for private(ndsArray,fn,gp)
	for(int p=0;p<nmpmsNR;p++)
	{	MPMBase *mpmptr = mpm[p];
		try
		{	// get shape functions
			const ElementBase *elemRef = theElements[mpmptr->ElemID()];
			int *nds = ndsArray;
			elemRef->GetShapeFunctions(fn,&nds,mpmptr);
			int numnds = nds[0];
			
			// Update particle position and velocity
			const MaterialBase *matRef=theMaterials[mpmptr->MatID()];
			int matfld=matRef->GetField();
			
			// Allow material to override global settings
			gp.m = m;
			gp.gridAlpha = gridAlpha;
			gp.particleAlpha = matRef->GetMaterialDamping(particleAlpha);

			// extrapolate nodal velocity from grid to particle S v+(k)
			ZeroVector(&gp.Svtilde);
			
			if(m<=0)
			{	// acceleration on the particle or S a (for FLIP and XPIC)
				ZeroVector(&gp.Sacc);
				
				// XPIC(k>1) needs separate velocity extrapolation
				if(m<-1) ZeroVector(&gp.Svlumped);
			}
			
			// only two possible transport tasks
			double rate[2],value[2];
			if(transportTasks!=NULL)
			{	rate[0] = rate[1] = value[0] = value[1] = 0.;
			}
			int task;
			TransportTask *nextTransport;
			
			// Loop over nodes
			for(int i=1;i<=numnds;i++)
			{	// increment velocity and acceleraton
				NodalPoint *ndptr = nd[nds[i]];
				short vfld = (short)mpmptr->vfld[i];

				// increment
				ndptr->IncrementDelvaTask5(vfld,matfld,fn[i],&gp);
				
				// increment transport rates
				nextTransport=transportTasks;
				task=0;
				while(nextTransport!=NULL)
				{	value[task] += nextTransport->IncrementValueExtrap(ndptr,fn[i],vfld,matfld);
#ifdef TRANSPORT_FMPM
					if(!nextTransport->IsUsingTransportXPIC())
						rate[task] += nextTransport->IncrementTransportRate(ndptr,fn[i],vfld,matfld);
#else
					rate[task] += nextTransport->IncrementTransportRate(ndptr,fn[i],vfld,matfld);
#endif
					nextTransport = nextTransport->GetNextTransportTask();
					task++;
				}
			}
			
			// Update velocity and position
			mpmptr->MoveParticle(&gp);
			
			// update transport values
			ResidualStrains res;
			res.dT = res.dC = 0.;
			double dTcond = 0.,dTad = 0.;
			nextTransport=transportTasks;
			task=0;
			while(nextTransport!=NULL)
			{	// mechanics would need to add to store returned value on particle (or get delta from rate?)
				if(nextTransport == conduction)
				{	res.dT = nextTransport->GetDeltaValue(mpmptr,value[task]);
#ifdef TRANSPORT_FMPM
					dTcond = nextTransport->IsUsingTransportXPIC() ? res.dT : rate[task]*timestep;
#else
					dTcond = rate[task]*timestep;
#endif
				}
				else
					res.dC = nextTransport->GetDeltaValue(mpmptr,value[task]);
				nextTransport=nextTransport->MoveTransportValue(mpmptr,timestep,rate[task],value[task]);
				task++;
			}
			
			// energy coupling here adds adiabatic temperature rise
			if(ConductionTask::adiabatic)
			{	dTad = mpmptr->GetClear_dTad();						// in K
				mpmptr->pTemperature += dTad;						// in K
				mpmptr->pPreviousTemperature += dTad;				// in K
				res.dT += dTad;
			}
			
			// for heat energy and entropy
			if(ConductionTask::active)
			{	double dq = matRef->GetHeatCapacity(mpmptr)*dTcond;
				mpmptr->AddHeatEnergy(dq);
				mpmptr->AddEntropy(dq,mpmptr->pPreviousTemperature);
			}
			else
			{	// when conduction off, update previous temp here
				res.dT = mpmptr->pTemperature - mpmptr->pPreviousTemperature;
				mpmptr->pPreviousTemperature = mpmptr->pTemperature;
			}
			
#ifdef POROELASTICITY
			// poroelasticity update
			if(fmobj->HasPoroelasticity())
			{	double dpud = mpmptr->GetClear_dpud();					// in MPa
				mpmptr->pConcentration += dpud;							// in MPa
				
				// update previous and residual change
				mpmptr->pPreviousConcentration += dpud;				// in MPa
				res.dC += dpud;
				
				if(mpmptr->pPreviousConcentration<0.)
				{	// do not let pPreviousConcentration go negative
					// if C-dC>0 (previous dC), set dC to -previous, othersize zero
					res.dC = -fmax(mpmptr->pPreviousConcentration-res.dC,0.);
					mpmptr->pPreviousConcentration = 0.;
				}
				else if(mpmptr->pPreviousConcentration-res.dC<0.)
				{	// adjust dC if was negative, but now position to be new positive value
					res.dC = mpmptr->pPreviousConcentration;
				}
			}
#endif
			
			// for generalized plane stress or strain, increment szz if needed
			if(fmobj->np==PLANE_STRESS_MPM || fmobj->np==PLANE_STRAIN_MPM)
			{	res.doopse = mpmptr->oopIncrement;
				mpmptr->oopIncrement = 0.;
			}
			
			// store increments on the particles
			mpmptr->dTrans = res;
			
		}
		catch(CommonException& err)
		{	if(upErr==NULL)
			{
#pragma omp critical (error)
				upErr = new CommonException(err);
			}
		}
		catch(CommonException* err)
		{	if(upErr==NULL)
			{
#pragma omp critical (error)
				upErr = new CommonException(*err);
			}
		}
		catch(std::bad_alloc&)
		{	if(upErr==NULL)
			{
#pragma omp critical (error)
				upErr = new CommonException("Memory error","UpdateParticlesTask::Execute");
			}
		}
		catch(...)
		{	if(upErr==NULL)
			{
#pragma omp critical (error)
				upErr = new CommonException("Unexpected error","UpdateParticlesTask::Execute");
			}
		}
	}
	
	// throw any errors
	if(upErr!=NULL) throw *upErr;
    
    // rigid materials move at their current velocity
    for(int p=nmpmsRB;p<nmpms;p++)
    {	mpm[p]->MovePosition(timestep);
    }
    
    return true;
}
