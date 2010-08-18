/******************************************************************************************
	GridForcesTask.cpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	Find all forces on the grid including internal forces (from particle stress)
		external forces, and body forces. Add equivalent forces for transport
		tasks.
 
	After main loop, get traction forces on cracks, add crack tip heating,
		and imperfect interface forces on cracks and track interface energy.
 
	Sum all forces only with external damping on the nodes
 
	Finally reconcile forces with boundary conditions. Do same for transport
		tasks.
 
	Input Variables
		mpm[]->ncpos, sp, pFext
		nd[]->pk
		bdyFrc.damping
 
	Output Variables
		theMaterials->LoadTransportProperties() - changed if depend on particle state
		mvf[]->fint, fext, ftot
		nd[]->fdiff, fcond
		mpm[]->dispEnergy
******************************************************************************************/

#include "NairnMPM_Class/GridForcesTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/TransportTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Cracks/CrackNode.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"

#pragma mark CONSTRUCTORS

GridForcesTask::GridForcesTask(const char *name) : MPMTask(name)
{
	// zero this in case in 2D analysis
	for(int i=0;i<MaxShapeNds;i++) zDeriv[i]=0.;
}

#pragma mark REQUIRED METHODS

// Get total grid point forces (except external forces)
void GridForcesTask::Execute(void)
{
#ifdef _PROFILE_TASKS_
	double beginTime=fmobj->CPUTime();
#endif

	int numnds,nds[MaxShapeNds];
	MaterialBase *matID;
	double xfrc,yfrc,fn[MaxShapeNds],xDeriv[MaxShapeNds],yDeriv[MaxShapeNds];
	TransportTask *nextTransport;
	Vector theFrc;
	NodalPoint *ndptr;
	short vfld;
	
	// loop over particles
    for(int p=0;p<nmpms;p++)
	{	MPMBase *mpmptr=mpm[p];							// material point pointer
		matID=theMaterials[mpmptr->MatID()];			// material class object
		
		// skip if material is rigid
		if(matID->Rigid()) continue;
		
		// get transport tensors (if needed)
		if(transportTasks!=NULL)
			matID->LoadTransportProps(mpmptr,fmobj->np);
		
        double mp=mpmptr->mp;					// in g
		int matfld=matID->GetField();		// material field
		
        // find shape functions and derviatives
 		theElements[mpmptr->ElemID()]->
				GetShapeGradients(&numnds,fn,nds,mpmptr->GetNcpos(),xDeriv,yDeriv,zDeriv);
		
        // Add particle property to each node in the element
        for(int i=1;i<=numnds;i++)
		{	vfld=(short)mpmptr->vfld[i];					// crack velocity field to use
			ndptr=nd[nds[i]];					// nodal point pointer
			
            // internal force vector (in g mm/sec^2 or micro N)
			//	(note: stress is specific stress in units N/m^2 cm^3/g
			//	Multiply by 1000 to make it mm/sec^2)
			mpmptr->Fint(theFrc,xDeriv[i],yDeriv[i],zDeriv[i]);
			ndptr->AddFintTask3(vfld,matfld,&theFrc);
            
            // body forces (not 3D yet)
			if(bodyFrc.GetGravity(&xfrc,&yfrc))
			{	theFrc=MakeVector(mp*fn[i]*xfrc,mp*fn[i]*yfrc,0.);
				ndptr->AddFintTask3(vfld,matfld,&theFrc);
			}
            
			// external force vector
			mpmptr->Fext(theFrc,fn[i]);
            ndptr->AddFextTask3(vfld,matfld,&theFrc);
			
			// transport forces
			nextTransport=transportTasks;
			while(nextTransport!=NULL)
				nextTransport=nextTransport->AddForces(ndptr,mpmptr,fn[i],xDeriv[i],yDeriv[i],zDeriv[i]);
        }
		
		// clear coupled dissipated energy
		if(ConductionTask::energyCoupling) mpmptr->SetDispEnergy(0.);
	}
	
	// traction law forces add to mvf[]->fext
	if(fmobj->hasTractionCracks)
	{	CrackHeader *nextCrack=firstCrack;
		while(nextCrack!=NULL)
		{	nextCrack->TractionFext();
			nextCrack=(CrackHeader *)nextCrack->GetNextObject();
		}
	}
	
	// crack tip heating adds to fcond
	if(conduction) conduction->AddCrackTipHeating();
	
	// interface forces added to mvf[]->fint and track interface energy
	CrackNode::InterfaceOnKnownNodes();
    
	// Find grid total forces with external damping
	double damping=bodyFrc.GetDamping();
    for(int i=1;i<=nnodes;i++)
		nd[i]->CalcFtotTask3(damping);
	
    // Imposed BCs on ftot to get correct grid BCs for velocity
    NodalVelBC::ConsistentGridForces();
	
	// Do similar to tranport property BCs
	nextTransport=transportTasks;
	while(nextTransport!=NULL)
		nextTransport=nextTransport->SetTransportForceBCs(timestep);
    
#ifdef _PROFILE_TASKS_
	totalTaskTime+=fmobj->CPUTime()-beginTime;
#endif
}
