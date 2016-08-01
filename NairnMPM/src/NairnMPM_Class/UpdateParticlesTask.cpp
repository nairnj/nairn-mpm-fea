/********************************************************************************
	UpdateParticlesTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	Using new grid results to update particle velocity, position, temperature
		and concentration. For rigid particles, however, just update postion
********************************************************************************/

#include "NairnMPM_Class/UpdateParticlesTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Custom_Tasks/TransportTask.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark CONSTRUCTORS

UpdateParticlesTask::UpdateParticlesTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Update particle position, velocity, temp, and conc
void UpdateParticlesTask::Execute(void)
{
	CommonException *upErr = NULL;
	int numnds,ndsArray[maxShapeNodes];
	double fn[maxShapeNodes];
    
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
#pragma omp parallel for private(numnds,ndsArray,fn)
    for(int p=0;p<nmpmsNR;p++)
	{	MPMBase *mpmptr = mpm[p];
		
		try
		{	// get shape functions
			const ElementBase *elemRef = theElements[mpmptr->ElemID()];
			int *nds = ndsArray;
			elemRef->GetShapeFunctions(fn,&nds,mpmptr);
			numnds = nds[0];
			
			// Update particle position and velocity
			const MaterialBase *matRef=theMaterials[mpmptr->MatID()];
			int matfld=matRef->GetField();
            
			// Allow material to override global settings
            double localParticleAlpha = particleAlpha;
            double localGridAlpha = gridAlpha;
            matRef->GetMaterialDamping(localParticleAlpha,localGridAlpha,nonPICGridAlpha,globalPIC);
			
			// data structure for extrapolations
			GridToParticleExtrap *gp = new GridToParticleExtrap;
			
			// acceleration on the particle
			gp->acc = mpmptr->GetAcc();
			ZeroVector(gp->acc);
			
			// extrapolate nodal velocity from grid to particle
			ZeroVector(&gp->vgpnp1);
			
			// only two possible transport tasks
			double rate[2];
			rate[0] = rate[1] = 0.;
			int task;
			TransportTask *nextTransport;
			
			// Loop over nodes
			for(int i=1;i<=numnds;i++)
			{	// increment velocity and acceleraton
				const NodalPoint *ndptr = nd[nds[i]];
                short vfld = (short)mpmptr->vfld[i];
				
				// increment
				ndptr->IncrementDelvaTask5(vfld,matfld,fn[i],gp);

#ifdef CHECK_NAN
				// conditionally compiled check for nan velocities
                if(gp->vgpnp1.x!=gp->vgpnp1.x || gp->vgpnp1.y!=gp->vgpnp1.y || gp->vgpnp1.z!=gp->vgpnp1.z)
                {
#pragma omp critical (output)
					{	cout << "\n# UpdateParticlesTask::Execute: bad material velocity field for vfld=" << vfld << "matfld=" << matfld << " fn[i]=" << fn[i] << endl;;
						PrintVector("#  Particle velocity vgpn1 = ",&gp->vgpnp1);
						cout << endl;
						ndptr->Describe();
					}
                }
#endif
				
				// increment transport rates
				nextTransport=transportTasks;
				task=0;
				while(nextTransport!=NULL)
					nextTransport=nextTransport->IncrementTransportRate(ndptr,fn[i],rate[task++]);
			}

			// Find grid damping acceleration parts =  ag*Vgp(n) = ag*(Vgp(n+1) - Agp(n)*dt)
			Vector accExtra = gp->vgpnp1;
			AddScaledVector(&accExtra, gp->acc, -timestep);
			ScaleVector(&accExtra,localGridAlpha);
			
			// update position, and must be before velocity update because updates need initial velocity
            // This section does second order update
			mpmptr->MovePosition(timestep,&gp->vgpnp1,&accExtra,localParticleAlpha);
			
			// update velocity in mm/sec
			mpmptr->MoveVelocity(timestep,&accExtra);
			
			// update transport values
			nextTransport=transportTasks;
			task=0;
			while(nextTransport!=NULL)
				nextTransport=nextTransport->MoveTransportValue(mpmptr,timestep,rate[task++]);
			
			// thermal ramp
			thermal.UpdateParticleTemperature(&mpmptr->pTemperature,timestep);

			// energy coupling here adds adiabtic temperature rise
			if(ConductionTask::adiabatic)
			{	double dTad = mpmptr->GetBufferClear_dTad();			// in K
				mpmptr->pTemperature += dTad;							// in K
			}
			
			// delete grid to particle extrap data
			delete gp;
			
		}
		catch(CommonException err)
		{	if(upErr==NULL)
			{
#pragma omp critical (error)
				upErr = new CommonException(err);
			}
		}
    }
	
	// throw any errors
	if(upErr!=NULL) throw *upErr;
    
    // rigid materials move at their current velocity
    for(int p=nmpmsNR;p<nmpms;p++)
    {	mpm[p]->MovePosition(timestep);
    }
}
