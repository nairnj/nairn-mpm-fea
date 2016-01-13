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
	
	int numnds,nds[maxShapeNodes];
	double fn[maxShapeNodes];
	Vector vgpnp1;
    
    // Damping terms on the grid or on the particles
    //      particleAlpha   =  alpha(PIC)/dt + pdamping(t)
    //      gridAlpha       = -alpha(PIC)/dt + damping(t)
    double particleAlpha = bodyFrc.GetParticleDamping(mtime);
	double gridAlpha = bodyFrc.GetDamping(mtime);

	//		nonPICGridAlpha = damping(t)
	//		globalPIC       = alpha(PIC)/dt
    double nonPICGridAlpha = bodyFrc.GetNonPICDamping(mtime);
    double globalPIC = bodyFrc.GetPICDamping();

    // Update particle position, velocity, temp, and conc
#pragma omp parallel for private(numnds,nds,fn,vgpnp1)
    for(int p=0;p<nmpmsNR;p++)
	{	MPMBase *mpmptr = mpm[p];
		
		try
		{	// get shape functions
			const ElementBase *elemRef = theElements[mpmptr->ElemID()];
			elemRef->GetShapeFunctions(&numnds,fn,nds,mpmptr);
			
			// Update particle position and velocity
			const MaterialBase *matRef=theMaterials[mpmptr->MatID()];
			int matfld=matRef->GetField();
            
			// Allow material to override global settings
            double localParticleAlpha = particleAlpha;
            double localGridAlpha = gridAlpha;
            matRef->GetMaterialDamping(localParticleAlpha,localGridAlpha,nonPICGridAlpha,globalPIC);
			
			Vector *acc=mpmptr->GetAcc();
			ZeroVector(acc);
			ZeroVector(&vgpnp1);
			double rate[2];         // only two possible transport tasks
			rate[0] = rate[1] = 0.;
			int task;
			TransportTask *nextTransport;
            short vfld;
			
			// Loop over nodes
			for(int i=1;i<=numnds;i++)
			{	// increment velocity and acceleraton
				const NodalPoint *ndptr = nd[nds[i]];
                vfld = (short)mpmptr->vfld[i];
				ndptr->IncrementDelvaTask5(vfld,matfld,fn[i],&vgpnp1,acc);

#ifdef CHECK_NAN
                if(vgpnp1.x!=vgpnp1.x || vgpnp1.y!=vgpnp1.y || vgpnp1.z!=vgpnp1.z)
                {
#pragma omp critical (output)
					{	cout << "\n# UpdateParticlesTask::Execute: bad material velocity field for vfld = " << vfld;
						PrintVector(" vgpn1 = ",&vgpnp1);
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
			
			// Find vgpn
			Vector vgpn = vgpnp1;
			AddScaledVector(&vgpn,acc,-timestep);
            
			// find effective grid acceleration and velocity
			AddScaledVector(acc,&vgpn,-localGridAlpha);
			AddScaledVector(&vgpnp1,&vgpn,-timestep*localGridAlpha);
			
			// update position, and must be before velocity update because updates need initial velocity
            // This section does second order update
			mpmptr->MovePosition(timestep,&vgpnp1,0.5*timestep,localParticleAlpha);

			// update velocity in mm/sec
			mpmptr->MoveVelocity(timestep,localParticleAlpha);
			
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
				mpmptr->pTemperature += dTad;                       // in K
			}
			
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
    {	mpm[p]->MovePosition(timestep,&mpm[p]->vel,0.,0.);
    }
}
