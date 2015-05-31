/********************************************************************************
	PostForcesTask.cpp
	nairn-mpm-fea

	Created by John Nairn on March 8, 2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.

	This task has various task that are done after extrapolating
	forces to the grid:

	1. Add traction BC nodal forces
	2. Add crack tractions to nodal forces
	3. Add crack tip heating
	4. Add imperfectinterface forces to nodes
	5. Add gravity and body forces
	6. Make nodes with velocity BCs have consistent forces
	7. Make transport tasks forces consistent with transport nodal BCs
********************************************************************************/

#include "NairnMPM_Class/PostForcesTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Nodes/MaterialInterfaceNode.hpp"
#include "Cracks/CrackNode.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "Boundary_Conditions/MatPtTractionBC.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"

#pragma mark CONSTRUCTORS

PostForcesTask::PostForcesTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get mass matrix, find dimensionless particle locations,
//	and find grid momenta
void PostForcesTask::Execute(void)
{
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
	
	// Add crack tip heating adds to conduction force
	if(conduction) conduction->AddCrackTipHeating();
	
	// Add interface forces to velocity fields and track total interface energy
    NodalPoint::interfaceEnergy=0.;
    CrackNode::CrackInterfaceOnKnownNodes();
    MaterialInterfaceNode::MaterialInterfaceOnKnownNodes();
    
	// Add gravity and body forces (if any are present)
	// Note: If ever need to implement body force that depend on particle state (stress, strain, etc.)
	//			then move the body force addition into GridForcesTask loop where gravity is commented out
    // When used to keep Fext, this section would all add fint and fext to get ftot (and it was always needed)
	Vector gridBodyForce;
	if(bodyFrc.gravity || bodyFrc.hasGridBodyForce)
	{   for(int i=1;i<=nnodes;i++)
		{	NodalPoint *ndptr = nd[i];
			bodyFrc.GetGridBodyForce(&gridBodyForce,ndptr,mtime);
			ndptr->AddGravityAndBodyForceTask3(&gridBodyForce);
		}
	}
	
    // Impose BCs on ftot to get correct grid BCs for velocity
    NodalVelBC::ConsistentGridForces();
	
	// Do similar to transport property BCs (not parallel because small and possible use of function/global variables)
	TransportTask *nextTransport=transportTasks;
	while(nextTransport!=NULL)
		nextTransport=nextTransport->SetTransportForceBCs(timestep);
}
