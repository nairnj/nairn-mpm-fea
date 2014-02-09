/********************************************************************************
	UpdateParticlesTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	Using new grid results to update particle velocity, position, temperature
		and concentration. For rigid particles, however, just update postion
 
	Input Variables
		mvf[]->pk, ftot, mass
		nd[]->fdiff, fcond
		bdyFrc.alpha
 
	Output Variabels
		mpm[]->pos, vel, acc, pConcentration, pTemperature
		bdyFrc.alpha, kineticEnergy, totalMass
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
	Vector delv;
	
    // Update particle position, velocity, temp, and conc
#pragma omp parallel for private(numnds,nds,fn,delv)
    for(int p=0;p<nmpmsNR;p++)
	{	MPMBase *mpmptr = mpm[p];
		
		try
		{	// get shape functions
			const ElementBase *elemRef = theElements[mpmptr->ElemID()];
			elemRef->GetShapeFunctions(&numnds,fn,nds,mpmptr);
			
			// Update particle position and velocity
			const MaterialBase *matRef=theMaterials[mpmptr->MatID()];
			int matfld=matRef->GetField();
			
			Vector *acc=mpmptr->GetAcc();
			ZeroVector(acc);
			ZeroVector(&delv);
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
				ndptr->IncrementDelvaTask5(vfld,matfld,fn[i],&delv,acc);

#ifdef CHECK_NAN
                if(delv.x!=delv.x || delv.y!=delv.y || delv.z!=delv.z)
                {   cout << "\n# UpdateParticlesTask::Execute: bad material velocity field for vfld = " << vfld << endl;
                    ndptr->Describe();
                }
#endif
				
				// increment transport rates
				nextTransport=transportTasks;
				task=0;
				while(nextTransport!=NULL)
					nextTransport=nextTransport->IncrementTransportRate(ndptr,fn[i],rate[task++]);
			}
			
			// update position in mm and velocity in mm/sec
			mpmptr->MovePosition(timestep,&delv);
			mpmptr->MoveVelocity(timestep,bodyFrc.GetAlpha(),&delv);
			
			// update transport values
			nextTransport=transportTasks;
			task=0;
			while(nextTransport!=NULL)
				nextTransport=nextTransport->MoveTransportValue(mpmptr,timestep,rate[task++]);
			
			// thermal ramp
			thermal.UpdateParticleTemperature(&mpmptr->pTemperature,timestep);
			
			// energy coupling here if conduction not doing it
			if(!ConductionTask::active)
			{	if(ConductionTask::adiabatic)
				{	double energy = mpmptr->GetDispEnergy();				// in uJ/g
					double Cv=1000.*matRef->GetHeatCapacity(mpmptr);		// in uJ/(g-K)
					mpmptr->pTemperature += energy/Cv;                      // in K
				}
				mpmptr->SetDispEnergy(0.);
			}
		}
		catch(CommonException err)
		{	if(upErr==NULL)
			{
#pragma omp critical
				upErr = new CommonException(err);
			}
		}
    }
	
	// throw any errors
	if(upErr!=NULL) throw *upErr;
    
    // rigid materials move at their current velocity
    for(int p=nmpmsNR;p<nmpms;p++)
    {	mpm[p]->MovePosition(timestep,&mpm[p]->vel);
    }
}
