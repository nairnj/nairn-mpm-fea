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

#pragma mark CONSTRUCTORS

UpdateParticlesTask::UpdateParticlesTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Update particle position, velocity, temp, and conc
// throws CommonException()
void UpdateParticlesTask::Execute(void)
{
	CommonException *upErr = NULL;
#ifdef CONST_ARRAYS
	int ndsArray[MAX_SHAPE_NODES];
	double fn[MAX_SHAPE_NODES];
#else
	int ndsArray[maxShapeNodes];
	double fn[maxShapeNodes];
#endif

    // Damping terms on the grid or on the particles
    //      particleAlpha   =  (1-beta)/dt + pdamping(t)
    //      gridAlpha       = -m*(1-beta)/dt + damping(t)
    double particleAlpha = bodyFrc.GetParticleDamping(mtime);
	double gridAlpha = bodyFrc.GetDamping(mtime);

	//		nonPICGridAlpha = damping(t)
	//		globalPIC       = (1-beta)/dt
    double nonPICGridAlpha = bodyFrc.GetNonPICDamping(mtime);
    double globalPIC = bodyFrc.GetPICDamping();

    // Update particle position, velocity, temp, and conc
#pragma omp parallel for private(ndsArray,fn)
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
            double localParticleAlpha = particleAlpha;
            double localGridAlpha = gridAlpha;
            matRef->GetMaterialDamping(localParticleAlpha,localGridAlpha,nonPICGridAlpha,globalPIC);
			
			// data structure for extrapolations
			GridToParticleExtrap *gp = new GridToParticleExtrap;
			
			// acceleration on the particle or S a
			gp->acc = mpmptr->GetAcc();
			ZeroVector(gp->acc);
			
			// extrapolate nodal velocity from grid to particle S v+
			ZeroVector(&gp->vgpnp1);
			
			// only two possible transport tasks
			double rate[2],value[2];
			rate[0] = rate[1] = value[0] = value[1] = 0.;
			int task;
			TransportTask *nextTransport;
			
			// Loop over nodes
			for(int i=1;i<=numnds;i++)
			{	// increment velocity and acceleraton
				NodalPoint *ndptr = nd[nds[i]];
                short vfld = (short)mpmptr->vfld[i];
				
				// increment
				ndptr->IncrementDelvaTask5(vfld,matfld,fn[i],gp);

				// increment transport rates
				nextTransport=transportTasks;
				task=0;
				while(nextTransport!=NULL)
				{	value[task] += nextTransport->IncrementValueExtrap(ndptr,fn[i],vfld,matfld);
					nextTransport = nextTransport->IncrementTransportRate(ndptr,fn[i],rate[task++],vfld,matfld);
				}
			}
			
			// Find initial velocity vgp = vgpnp1 - a dt (may be better to not need this calculation)
			Vector vgpn = gp->vgpnp1;
			AddScaledVector(&vgpn, gp->acc, -timestep);

			// get extra particle acceleration Aextra = -ag*Vgp(n) - ap*Vp - m(1-beta)Sv^*/dt
			// (the -ap*Vp term is done in MovePosition, others done here)
			Vector accExtra = SetScaledVector(&vgpn,-localGridAlpha);
			
			// update position, and must be before velocity update because updates need initial velocity
            // This section does second order update
			mpmptr->MovePosition(timestep,&gp->vgpnp1,&accExtra,localParticleAlpha);
			
			// update velocity in mm/sec
			mpmptr->MoveVelocity(timestep);
			
			// update transport values
			ResidualStrains res;
			res.dT = res.dC = 0.;
			double dTcond = 0.,dTad = 0.;
			nextTransport=transportTasks;
			task=0;
			while(nextTransport!=NULL)
			{	// mechanics would need to add to store returned value on particle (or get delta from rate?)
				if(nextTransport == conduction)
				{	dTcond = rate[task]*timestep;
					res.dT = nextTransport->GetDeltaValue(mpmptr,value[task]);;
				}
				else
					res.dC = nextTransport->GetDeltaValue(mpmptr,value[task]);;
				nextTransport=nextTransport->MoveTransportValue(mpmptr,timestep,rate[task++]);
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
			
			// for generalized plane stress or strain, increment szz if needed
			if(fmobj->np==PLANE_STRESS_MPM || fmobj->np==PLANE_STRAIN_MPM)
			{	res.doopse = mpmptr->oopIncrement;
				mpmptr->oopIncrement = 0.;
			}
			
			// store increments on the particles
			mpmptr->dTrans = res;
			
			// delete grid to particle extrap data
			delete gp;
			
		}
		catch(CommonException& err)
		{	if(upErr==NULL)
			{
#pragma omp critical (error)
				upErr = new CommonException(err);
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
}
