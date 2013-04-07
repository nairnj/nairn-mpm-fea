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
******************************************************************************************/

#include "NairnMPM_Class/GridForcesTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/TransportTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Nodes/MaterialInterfaceNode.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Cracks/CrackNode.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Boundary_Conditions/MatPtTractionBC.hpp"
#ifdef LOG_PROGRESS
#include "System/ArchiveData.hpp"
#endif


#pragma mark CONSTRUCTORS

GridForcesTask::GridForcesTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get total grid point forces (except external forces)
void GridForcesTask::Execute(void)
{
	// need to be private in threads
	TransportProperties t;
	int numnds,nds[maxShapeNodes];
	double fn[maxShapeNodes],xDeriv[maxShapeNodes],yDeriv[maxShapeNodes],zDeriv[maxShapeNodes];
	
	// loop over particles (if made parallel, very close due to needed critical or atomic sections
#pragma omp parallel for private(t,numnds,nds,fn,xDeriv,yDeriv,zDeriv)
    for(int p=0;p<nmpmsNR;p++)
	{	MPMBase *mpmptr=mpm[p];											// material point pointer
		const MaterialBase *matref = theMaterials[mpmptr->MatID()];		// material class (read only)
		
		// skip if material is rigid (and comes before last nonrigid one)
		if(matref->Rigid()) continue;
		
		// get transport tensors (if needed)
		if(transportTasks!=NULL)
			matref->GetTransportProps(mpmptr,fmobj->np,&t);
		
        double mp = mpmptr->mp;					// in g
		//int matfld = matref->GetField();           // material field
		
        // find shape functions and derviatives
		const ElementBase *elemref = theElements[mpmptr->ElemID()];
 		elemref->GetShapeGradients(&numnds,fn,nds,mpmptr->GetNcpos(),xDeriv,yDeriv,zDeriv,mpmptr);
		
        // Add particle property to each node in the element
		mpmptr->vfld[0] = numnds;
        for(int i=1;i<=numnds;i++)
		{	//short vfld = (short)mpmptr->vfld[i];						// crack velocity field to use
			//NodalPoint *ndptr = nd[nds[i]];							// nodal point pointer
			
            // total force vector (in g mm/sec^2 or micro N)
			//	(note: stress is specific stress in units N/m^2 cm^3/g
			//	Multiply by 1000 to make it mm/sec^2)
			Vector theFrc;
			mpmptr->GetFint(theFrc,xDeriv[i],yDeriv[i],zDeriv[i]);
            
            // add body forces
			bodyFrc.AddGravity(mp,fn[i],&theFrc);
            
			// get external force vector and add to velocity field
			mpmptr->AddFext(theFrc,fn[i]);
			
			// load into force buffer
			int i0 = i-1;
			mpmptr->gFrc[i0].nodeNum = nds[i];
			mpmptr->gFrc[i0].forces[0] = theFrc.x;
			mpmptr->gFrc[i0].forces[1] = theFrc.y;
			mpmptr->gFrc[i0].forces[2] = theFrc.z;
			
			// this is critical section that will change nodal values
			
			/*
			// Now add total force to the node
			ndptr->AddFtotTask3(vfld,matfld,&theFrc);
		
			// transport forces
			TransportTask *nextTransport=transportTasks;
			while(nextTransport!=NULL)
				nextTransport=nextTransport->AddForces(ndptr,mpmptr,fn[i],xDeriv[i],yDeriv[i],zDeriv[i],&t);
			*/
        }
		
		// clear coupled dissipated energy if in conduction becaouse done with it this time step
		if(ConductionTask::active) mpmptr->SetDispEnergy(0.);
	}
	
	for(int p=0;p<nmpmsNR;p++)
	{	MPMBase *mpmptr=mpm[p];											// material point pointer
		const MaterialBase *matref = theMaterials[mpmptr->MatID()];		// material class (read only)
		
		// skip if material is rigid (and comes before last nonrigid one)
		if(matref->Rigid()) continue;
		int matfld = matref->GetField();           // material field
		
        // Add particle property to each node in the element
        for(int i=0;i<(int)mpmptr->vfld[0];i++)
		{	short vfld = (short)mpmptr->vfld[i+1];						// crack velocity field to use
			NodalPoint *ndptr = nd[mpmptr->gFrc[i].nodeNum];					// nodal point pointer
			Vector theFrc = MakeVector(mpmptr->gFrc[i].forces[0],mpmptr->gFrc[i].forces[1],mpmptr->gFrc[i].forces[2]);
			ndptr->AddFtotTask3(vfld,matfld,&theFrc);
		}
	}
	
	// Add traction BCs on particles
	MatPtTractionBC::SetParticleSurfaceTractions(mtime);
	
	// Add traction law forces to velocity fields
	if(fmobj->hasTractionCracks)
	{	CrackHeader *nextCrack=firstCrack;
		while(nextCrack!=NULL)
		{	nextCrack->AddTractionForce();
			nextCrack=(CrackHeader *)nextCrack->GetNextObject();
		}
	}
	
	// Add crack tip heating adds to fcond or conduction force
	if(conduction) conduction->AddCrackTipHeating();
	
	// Add interface forces added to velocity fields and track total interface energy
    NodalPoint::interfaceEnergy=0.;
    CrackNode::InterfaceOnKnownNodes();
    MaterialInterfaceNode::InterfaceOnKnownNodes();
    
	// Add grid damping forces
	if(bodyFrc.useDamping)
	{	double damping = bodyFrc.GetDamping(mtime);		// could move inside loop and make function of nodal position too
    	for(int i=1;i<=nnodes;i++)
			nd[i]->AddGridDampingTask3(damping);
	}
	
    // Imposed BCs on ftot to get correct grid BCs for velocity
    NodalVelBC::ConsistentGridForces();
	
	// Do similar to transport property BCs
	TransportTask *nextTransport=transportTasks;
	while(nextTransport!=NULL)
		nextTransport=nextTransport->SetTransportForceBCs(timestep);
    
}
