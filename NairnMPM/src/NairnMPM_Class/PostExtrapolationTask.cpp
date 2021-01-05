/********************************************************************************
	PostExtrapolationTask.cpp
	nairn-mpm-fea

	Created by John Nairn on March 6, 2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.
 
	This task has various task that are done after extrapolating
	mass and momentum to the grid. The tasks are:
	---------------------------------------------
	* Mirror crack fields (if needed)
	* Get total mass and count particles and copy no BC momentum
	* Find nodes with material or crack contact
	* Transport tasks get value
	  (also do for CVF and MVF if contact flow activated)
	* After main loop
 		- Find active nodes
 		- locate mirrored BCs
 		- material and crack contact
 		- impose velocity BCs
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
#include "Nodes/MaterialContactNode.hpp"
#include "NairnMPM_Class/UpdateMomentaTask.hpp"

#pragma mark CONSTRUCTORS

PostExtrapolationTask::PostExtrapolationTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get mass matrix, find dimensionless particle locations,
//	and find grid momenta
// throws CommonException()
bool PostExtrapolationTask::Execute(int taskOption)
{
	CommonException *massErr = NULL;

    // flag if combining particles are present (i.e., has cracks, multimaterial mode, and has rigid contact particles)
	bool mirrorIgnored = firstCrack!=NULL && fmobj->hasNoncrackingParticles;
	
	// First node pass does some calculations and gets transport values
	// If needed, find material contact nodes (if numberMaterials>1)
#pragma omp parallel
	{
		// variables for each thread
		MaterialContactNode *lastMCNode=NULL;
		CrackNode *lastCrackNode=NULL;
		
		// Each pass in this loop should be independent
#pragma omp for nowait
		for(int i=1;i<=nnodes;i++)
		{	// node reference
			NodalPoint *ndptr = nd[i];
			
			try
			{
				// mirror velocity fields for materials that ignore cracks
				if(mirrorIgnored)
					ndptr->MirrorIgnoredCrackFields();
				
				// Get total nodal masses and count materials if multimaterial mode
				if(ndptr->CalcTotalMassAndCount())
				{	// save multimaterial nodes that might have contact
					lastMCNode = new MaterialContactNode(ndptr,lastMCNode);
					ndptr->contactData = lastMCNode;
				}
				
				// if needed create crack contact node and added to this patches linked list
				if(firstCrack!=NULL)
				{	int hasFlags = ndptr->HasCrackContact();
					if(hasFlags)
					{	lastCrackNode = new CrackNode(ndptr,hasFlags,lastCrackNode);
					}
				}

				// get transport values on nodes
				TransportTask::GetTransportValues(ndptr);
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
	
		// as each thread exits above loop, insert into list (if was needed)
		if(lastMCNode != NULL)
		{
#pragma omp critical (linknodes)
			{	// add to node vector
				MaterialContactNode *prevMCNode = lastMCNode;
				while(prevMCNode!=NULL)
				{	MaterialContactNode::materialContactNodes.push_back(prevMCNode);
					prevMCNode = prevMCNode->GetPrevNode();
				}
			}
		}
		
		// as each thread exits above loop, insert into list (if was needed)
		if(lastCrackNode != NULL)
		{
#pragma omp critical (linknodes)
			{	// add to node vector
				CrackNode *prevCNode = lastCrackNode;
				while(prevCNode!=NULL)
				{	CrackNode::crackContactNodes.push_back(prevCNode);
					prevCNode = prevCNode->GetPrevNode();
				}
			}
		}
	}
	
	// throw any errors
	if(massErr!=NULL) throw *massErr;	// Post mass and momentum extrapolation calculations on nodes

	// Create list of active nodes
	int numActive = 0;
	for(int i=1;i<=nnodes;i++)
	{	if(nd[i]->NodeHasParticles())
			nda[++numActive] = i;
	}
	nda[0] = numActive;

	// locate BCs with reflected nodes
	if(firstRigidVelocityBC!=NULL)
	{   NodalVelBC *nextBC=firstRigidVelocityBC;
		while(nextBC!=NULL)
			nextBC = nextBC->SetMirroredVelBC(mtime);
	}
	
	// precalculate velocity BC values
	NodalVelBC::GridVelocityBCValues();
	
	// contact and grid velocity conditions
	UpdateMomentaTask::ContactAndMomentaBCs(MASS_MOMENTUM_CALL);

	// Impose transport BCs and extrapolate gradients to the particles
	TransportTask::TransportBCsAndGradients(mtime);
    
    return true;
}
