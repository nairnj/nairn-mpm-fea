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

#pragma mark CONSTRUCTORS

UpdateParticlesTask::UpdateParticlesTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Update particle position, velocity, temp, and conc
void UpdateParticlesTask::Execute(void)
{
#ifdef _PROFILE_TASKS_
	double beginTime=fmobj->CPUTime();
#endif

	int i,p,iel,numnds,matfld,nds[MaxShapeNds];
	MaterialBase *matID;
	TransportTask *nextTransport;
	double fn[MaxShapeNds];
	
    // Update particle position, velocity, temp, and conc
	Vector delv,*acc;
	bodyFrc.TrackAlpha();       // zero total kinectic energy
    for(p=0;p<nmpms;p++)
	{	matID=theMaterials[mpm[p]->MatID()];
		if(!matID->Rigid())
		{	// get shape functions
			iel=mpm[p]->ElemID();
			theElements[iel]->GetShapeFunctions(&numnds,fn,nds,mpm[p]->GetNcpos());
			
			// Update particle position and velocity
			matfld=matID->GetField();
			acc=mpm[p]->GetAcc();
			ZeroVector(acc);
			ZeroVector(&delv);
			nextTransport=transportTasks;
			while(nextTransport!=NULL)
				nextTransport=nextTransport->ZeroTransportRate();
			//int mpmvfld=mpm[p]->vfld[1];
			//bool sameField=TRUE;
			for(i=1;i<=numnds;i++)
			{	nd[nds[i]]->IncrementDelvaTask5((short)mpm[p]->vfld[i],matfld,fn[i],&delv,acc);
				//if(mpm[p]->vfld[i]!=mpmvfld) sameField=FALSE;
				nextTransport=transportTasks;
				while(nextTransport!=NULL)
					nextTransport=nextTransport->IncrementTransportRate(nd[nds[i]],fn[i]);
			}
			
			//if(!sameField) cout << "see " << p+1 << endl;
			
			// update position in mm and velocity in mm/sec
			mpm[p]->MovePosition(timestep,&delv);
			mpm[p]->MoveVelocity(timestep,bodyFrc.GetAlpha(),&delv);
			
			// update transport values
			nextTransport=transportTasks;
			while(nextTransport!=NULL)
				nextTransport=nextTransport->MoveTransportValue(mpm[p],timestep);
			
			// thermal ramp
			thermal.UpdateParticleTemperature(&mpm[p]->pTemperature,timestep);
			
			// update feedback coefficient
			bodyFrc.TrackAlpha(mpm[p]);
		}
		
		else
		{	// rigid materials at constant velocity
			mpm[p]->MovePosition(timestep,&mpm[p]->vel);
		}
	}
	
	// update damping coefficient
	bodyFrc.UpdateAlpha(timestep,mtime);
	
#ifdef _PROFILE_TASKS_
	totalTaskTime+=fmobj->CPUTime()-beginTime;
#endif
}
