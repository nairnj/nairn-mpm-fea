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
	
	// loop over non-rigid particles - this parallel part changes only particle p
	// forces are stored on buffer, which are sent to nodes in next non-parallel loop
    for(int p=0;p<nmpmsNR;p++)
	{	MPMBase *mpmptr=mpm[p];											// material point pointer
		const MaterialBase *matref = theMaterials[mpmptr->MatID()];		// material class (read only)
        int matfld = matref->GetField(); 
		
		// get transport tensors (if needed)
		if(transportTasks!=NULL)
			matref->GetTransportProps(mpmptr,fmobj->np,&t);
		
        // find shape functions and derviatives
		const ElementBase *elemref = theElements[mpmptr->ElemID()];
 		elemref->GetShapeGradients(&numnds,fn,nds,mpmptr->GetNcpos(),xDeriv,yDeriv,zDeriv,mpmptr);
		mpmptr->vfld[0] = numnds;			// save number of nodes
		
        // Add particle property to buffer on the material point (needed to allow parallel code)
        short vfld;
        NodalPoint *ndptr;
        for(int i=1;i<=numnds;i++)
		{	vfld = (short)mpmptr->vfld[i];					// crack velocity field to use
			ndptr = nd[nds[i]];								// nodal point pointer
            
			// total force vector = internal + external forces
			//	(in g mm/sec^2 or micro N)
            Vector theFrc;
			mpmptr->GetFintPlusFext(&theFrc,fn[i],xDeriv[i],yDeriv[i],zDeriv[i]);
            
            // add body forces
			bodyFrc.AddGravity(&theFrc,mpmptr->mp,fn[i]);
            
            // add to total force
            ndptr->AddFtotTask3(vfld,matfld,&theFrc);
			
			// transport forces
            TransportTask *nextTransport=transportTasks;
            while(nextTransport!=NULL)
                nextTransport=nextTransport->AddForces(ndptr,mpmptr,fn[i],xDeriv[i],yDeriv[i],zDeriv[i],&t);
        }
		
		// clear coupled dissipated energy if in conduction becaouse done with it this time step
		if(ConductionTask::active) mpmptr->SetDispEnergy(0.);
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
