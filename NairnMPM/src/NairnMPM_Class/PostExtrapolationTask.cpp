/********************************************************************************
	PostExtrapolationTask.cpp
	nairn-mpm-fea

	Created by John Nairn on March 6, 2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.
 
	This task has various task that are done after extrapolating
	mass and momentum to the grid. The tasks are:
	---------------------------------------------
	* Mirror crack fields (if needed)
	* Get total mass and count particles
	* Multimaterial contact (save the nodes for imperfect interfaces)
	* Crack contact (save the crack nodes)
	* Transport tasks get value
	  (also do for CVF and MVF if contact flow activated)
	* After main loop
		- Extrapolate temperature and concentration gradients to the grid
		- Set mirrored BCs
		- Copy no BC velocity and then impose all boundary conditions
		  (GridMomentumConditions())
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/PostExtrapolationTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Custom_Tasks/TransportTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "Cracks/CrackNode.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"

#pragma mark CONSTRUCTORS

PostExtrapolationTask::PostExtrapolationTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get mass matrix, find dimensionless particle locations,
//	and find grid momenta
// throws CommonException()
void PostExtrapolationTask::Execute(void)
{
	CommonException *massErr = NULL;
	
	bool combineRigid = firstCrack!=NULL && fmobj->multiMaterialMode && fmobj->hasRigidContactParticles;
	
	// Post mass and momentum extrapolation calculations on nodes
#pragma omp parallel
	{
		// variables for each thread
		CrackNode *firstCrackNode=NULL,*lastCrackNode=NULL;
		
		// Each pass in this loop should be independent
#pragma omp for nowait
		for(int i=1;i<=nnodes;i++)
		{	// node reference
			NodalPoint *ndptr = nd[i];
			
			try
            {	// combine rigid fields if necessary
                if(combineRigid)
                    ndptr->CopyRigidParticleField();
				// Get total nodal masses and count materials if multimaterial mode
				ndptr->CalcTotalMassAndCount();
				
				// multimaterial contact and interfaces
				if(fmobj->multiMaterialMode)
					ndptr->MaterialContactOnNode(timestep,MASS_MOMENTUM_CALL);
				
				// crack contact and interfaces
				if(firstCrack!=NULL)
					ndptr->CrackContact(MASS_MOMENTUM_CALL,0.,&firstCrackNode,&lastCrackNode);
				
				// get transport values on nodes
				TransportTask *nextTransport=transportTasks;
				while (nextTransport != NULL)
					nextTransport = nextTransport->GetTransportNodalValue(ndptr);
			}
			catch(std::bad_alloc&)
			{	if(massErr==NULL)
				{
#pragma omp critical (error)
					massErr = new CommonException("Memory error","PostExtrapolationTask::Execute");
				}
			}
			catch(...)
			{	if(massErr==NULL)
				{
#pragma omp critical (error)
					massErr = new CommonException("Unexpected error","PostExtrapolationTask::Execute");
				}
			}
		}
		
#pragma omp critical (linknodes)
		{
			// link up crack nodes
			if(lastCrackNode != NULL)
			{	if(CrackNode::currentCNode != NULL)
					firstCrackNode->SetPrevNode(CrackNode::currentCNode);
				CrackNode::currentCNode = lastCrackNode;
			}
		}
	}
	
	// throw any errors
	if(massErr!=NULL) throw *massErr;
	
#pragma mark ... IMPOSE BOUNDARY CONDITIONS
	
	// Impose transport BCs and extrapolate gradients to the particles
	TransportTask *nextTransport=transportTasks;
	while(nextTransport!=NULL)
	{   nextTransport->ImposeValueBCs(mtime);
		nextTransport = nextTransport->GetGradients(mtime);
	}
	
	// locate BCs with reflected nodes
	if(firstRigidVelocityBC!=NULL)
	{   NodalVelBC *nextBC=firstRigidVelocityBC;
		while(nextBC!=NULL)
			nextBC = nextBC->SetMirroredVelBC(mtime);
	}
	
	// used to call class methods for material contact and crack contact here
	// Impose velocity BCs
	NodalVelBC::GridMomentumConditions();
}
