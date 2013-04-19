/********************************************************************************
	UpdateParticlesTask.cpp
	NairnMPM

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

#pragma mark CONSTRUCTORS

UpdateParticlesTask::UpdateParticlesTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Update particle position, velocity, temp, and conc
void UpdateParticlesTask::Execute(void)
{
	int numnds,nds[maxShapeNodes];
	double fn[maxShapeNodes];
	Vector delv;
	
    // Update particle position, velocity, temp, and conc
#pragma omp parallel for private(numnds,nds,fn,delv)
    for(int p=0;p<nmpmsNR;p++)
	{	// get shape functions
        const ElementBase *elemRef = theElements[mpm[p]->ElemID()];
        elemRef->GetShapeFunctions(&numnds,fn,nds,mpm[p]->GetNcpos(),mpm[p]);
        
        // Update particle position and velocity
        const MaterialBase *matRef=theMaterials[mpm[p]->MatID()];
        int matfld=matRef->GetField();
        
        Vector *acc=mpm[p]->GetAcc();
        ZeroVector(acc);
        ZeroVector(&delv);
        double rate[2];         // only two possible transport tasks
        rate[0] = rate[1] = 0.;
        int task;
        TransportTask *nextTransport;
        
        // Loop over nodes
        for(int i=1;i<=numnds;i++)
        {	// increment velocity and acceleraton
            const NodalPoint *ndptr = nd[nds[i]];
            ndptr->IncrementDelvaTask5((short)mpm[p]->vfld[i],matfld,fn[i],&delv,acc);
            
            // increment transport rates
            nextTransport=transportTasks;
            task=0;
            while(nextTransport!=NULL)
                nextTransport=nextTransport->IncrementTransportRate(nd[nds[i]],fn[i],rate[task++]);
        }
        
        // update position in mm and velocity in mm/sec
        mpm[p]->MovePosition(timestep,&delv);
        mpm[p]->MoveVelocity(timestep,bodyFrc.GetAlpha(),&delv);
        
        // update transport values
        nextTransport=transportTasks;
        task=0;
        while(nextTransport!=NULL)
            nextTransport=nextTransport->MoveTransportValue(mpm[p],timestep,rate[task++]);
        
        // thermal ramp
        thermal.UpdateParticleTemperature(&mpm[p]->pTemperature,timestep);
        
        // energy coupling here if conduction not doing it
        if(!ConductionTask::active)
        {	if(ConductionTask::energyCoupling)
            {	double energy = mpm[p]->GetDispEnergy();									// in uJ/g
                double Cp=1000.*matRef->GetHeatCapacity(mpm[p]);		// in uJ/(g-K)
                mpm[p]->pTemperature += energy/Cp;			// in K
            }
            mpm[p]->SetDispEnergy(0.);
        }
    }
    
    // rigid materials move at their current velocity
    for(int p=nmpmsNR;p<nmpms;p++)
    {	mpm[p]->MovePosition(timestep,&mpm[p]->vel);
    }
}
