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
	* Update particle temperature, concentration, and other diffusion values
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
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "Global_Quantities/BodyForce.hpp"
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
	double rate[MAX_TRANSPORT_TASKS],value[MAX_TRANSPORT_TASKS],pic[MAX_TRANSPORT_TASKS];

	// data structure for extrapolations
	GridToParticleExtrap gp;
	
	// Damping terms on the grid or on the particles
	//      particleAlpha   =  pdamping(t)
	//      gridAlpha       = damping(t)
	double particleAlpha = bodyFrc.GetParticleDamping(mtime);
	double gridAlpha = bodyFrc.GetGridDamping(mtime);
	
	// Calculate velocities on the grid
	// For FLIP, FMPM(1), and XPIC(1) v_i will be v_i^{L+}
	// For FMPM(k>1) and XPIC(k>1) v(i) will be v(k) = m_k^{-1}p^+
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
	// m>0 is FMPM(m) while m<0 will be XPIC(-m) (m=0 is FLIP)
	if(!bodyFrc.UsingFMPM()) m = -m;
	
	// get grid transport values
	int xpicOrder = 0;
	if(XPICTransportTask!=NULL)
	{	xpicOrder = XPICTransportTask->GetXPICOrder();
		if(xpicOrder>1)
		{
#pragma omp parallel for
			for(int i=1;i<=*nda;i++)
				XPICTransportTask->CopyXStar(nd[nda[i]]);
			TransportTask::TransportGridBCs(mtime,timestep,UPDATE_GRID_STRAINS_CALL);
		}
	}

	// Update particle position, velocity, temp, conc, and other diffusion values
#pragma omp parallel for private(ndsArray,fn,gp,rate,value,pic)
	for(int p=0;p<nmpmsNR;p++)
	{	MPMBase *mpmptr = mpm[p];
		if(mpmptr->InReservoir()) continue;
		
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

			// extrapolate nodal velocity from grid to particle Sv
            // Sv+(lumped) FLIP/PIC, Sv+(k) FMPM, S(v(k)+a*dt) XPIC
			ZeroVector(&gp.Svk);
			
			if(m<=0)
			{	// acceleration on the particle or S a (for FLIP and XPIC)
				ZeroVector(&gp.Sacc);
			}
						
			// zero for each transport tasks
			for(int i=0;i<numTransport;i++)
			{	rate[i]=0.;
				value[i]=0.;
				pic[i]=0.;
			}
			
			// Loop over nodes for this particle and extrapolate grid results
			// to this particle
			int task;
			TransportTask *nextTransport;
			for(int i=1;i<=numnds;i++)
			{	// increment velocity and acceleraton
				NodalPoint *ndptr = nd[nds[i]];
				short vfld = (short)mpmptr->vfld[i];

				// increment
				ndptr->IncrementDelvaTask5(vfld,matfld,fn[i],&gp);
				
				// extraplate to the particle
				// FLIP: value and rate
				// FMPM(k): value, FLIP/FMPM(1): value and rate, FLIP/FMPM(k)>1; value, rate, and pic
				nextTransport=transportTasks;
				task=0;
				while(nextTransport!=NULL)
				{	// always get value
					value[task] += nextTransport->IncrementValueExtrap(ndptr,fn[i],vfld,matfld);
					double fractFMPM;
					if(nextTransport->IsUsingTransportXPIC(fractFMPM))
					{	// step includes FMPM(k), but FLIP only in fractMPM<1
						if(fractFMPM<1)
						{	// step blended FLIP/FMPM
							rate[task] += nextTransport->IncrementTransportRate(ndptr,fn[i],vfld,matfld);
							if(xpicOrder>1)
							{	// blended FLIP/FMPM(k>1)
								pic[task] += nextTransport->IncrementLumpedValueExtrap(ndptr,fn[i],vfld,matfld);
							}
						}
					}
					else
					{	// FLIP only
						rate[task] += nextTransport->IncrementTransportRate(ndptr,fn[i],vfld,matfld);
					}
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
				{	// Store grid extrapolation on the particle and get change from grid values
					double fractFMPM,gridValue=value[task];
					if(nextTransport->ShouldBlendFromGrid(fractFMPM))
						gridValue = fractFMPM*value[task] + (1.-fractFMPM)*pic[task];
					nextTransport->GetDeltaValue(mpmptr,gridValue,&res.dT);
					
					// dTcond is change in particle temperature due to conduction alone,
					//   which uses smoothed increment found above (changed to smoothed 10/25/2022)
					// It is used below to get heat energy, entropy, and phase transitions.
					dTcond = res.dT;
				}
				else
				{	// When coupled with mechanics, add particle source terms from
					// constitutive laws now
					int transNumber = ((DiffusionTask *)nextTransport)->GetNumber();
                    double dvalue;
                    if(mpmptr->GetClearParticleDiffusionSource(transNumber,dvalue))
                    {   value[task] += dvalue;
                        pic[task] += dvalue;
                        rate[task] += dvalue/timestep;
                    }

					// give task a chance to change rate and value(s)
					nextTransport->AdjustRateAndValue(mpmptr,value[task],rate[task],pic[task],timestep);
					
					// Store grid extrapolation on the partle and get change from grid values
					// Note that only base diffusion sets res.dC
					double fractFMPM,gridValue=value[task];
					if(nextTransport->ShouldBlendFromGrid(fractFMPM))
						gridValue = fractFMPM*value[task] + (1.-fractFMPM)*pic[task];
					nextTransport->GetDeltaValue(mpmptr,gridValue,&res.dC);
				}

				// change particle value, then on to next task
				nextTransport=nextTransport->MoveTransportValue(mpmptr,timestep,rate[task],value[task]);
				task++;
			}
			
			// energy coupling here adds adiabatic temperature rise
			if(ConductionTask::adiabatic)
			{	// increment temperatures and res.dT by Tad, but not dTcond
				dTad = mpmptr->GetClear_dTad();						// in K
				mpmptr->pTemperature += dTad;						// in K
				mpmptr->pPreviousTemperature += dTad;				// in K
				res.dT += dTad;
			}
			
			// for heat energy and entropy due to dTcond alone
			if(ConductionTask::active)
			{	// heat energy from dTcond only
				double cv = matRef->GetHeatCapacity(mpmptr);
				mpmptr->AddHeatEnergy(cv*dTcond);
				
				// Entropy for change due to dTcond alone (i.e. reversible)
				// Inreverible is handled when material calls IncrementHeatEnergy
				mpmptr->AddEntropy(cv,mpmptr->pPreviousTemperature-dTcond,mpmptr->pPreviousTemperature);
				
				// Isothermal mode subtracts dTad from dTcond, but that is handled
				// when material calls IncrementHeatEnergy
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
	{	if(mpm[p]->InReservoir()) continue;
		mpm[p]->MovePosition(timestep);
	}
    
    return true;
}
