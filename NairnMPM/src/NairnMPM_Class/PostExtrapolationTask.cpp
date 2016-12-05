/********************************************************************************
	PostExtrapolationTask.cpp
	nairn-mpm-fea

	Created by John Nairn on March 6, 2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.
 
	This task has various task that are done after extrapolating
	mass and momentum to the grid:
 
	1. Mirror crack fields
	3. Get total mass and count
	3. Multimaterial contact (save the nodes)
	4. Crack contact (save the nodes)
	5. transport tasks
	6. Impost all boundary conditions
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/PostExtrapolationTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Custom_Tasks/TransportTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "Cracks/CrackNode.hpp"
#include "Nodes/MaterialInterfaceNode.hpp"
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
		MaterialInterfaceNode *firstInterfaceNode=NULL,*lastInterfaceNode=NULL;
		
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
                
				// multimaterial contact
				if(fmobj->multiMaterialMode)
					ndptr->MaterialContactOnNode(timestep,MASS_MOMENTUM_CALL,&firstInterfaceNode,&lastInterfaceNode);
				
				// crack contact
				if(firstCrack!=NULL)
					ndptr->CrackContact(FALSE,0.,&firstCrackNode,&lastCrackNode);
				
				// get transport values on nodes
				TransportTask *nextTransport=transportTasks;
				while(nextTransport!=NULL)
					nextTransport = nextTransport->GetNodalValue(ndptr);
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
				firstCrackNode->SetPrevBC(CrackNode::currentCNode);
				CrackNode::currentCNode = lastCrackNode;
			}
			
			// link up interface nodes
			if(lastInterfaceNode != NULL)
			{	if(MaterialInterfaceNode::currentIntNode != NULL)
				firstInterfaceNode->SetPrevBC(MaterialInterfaceNode::currentIntNode);
				MaterialInterfaceNode::currentIntNode = lastInterfaceNode;
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
        //cout << "# Find Reflected Nodes" << endl;
        while(nextBC!=NULL)
            nextBC = nextBC->SetMirroredVelBC(mtime);
    }
	
	// used to call class methods for material contact and crack contact here
	// Impose velocity BCs
	NodalVelBC::GridMomentumConditions(TRUE);
}
