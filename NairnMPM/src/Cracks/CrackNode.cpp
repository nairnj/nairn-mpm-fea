/********************************************************************************
    CrackNode.cpp
    nairn-mpm-fea
    
    Created by John Nairn on 2/24/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	This object handles crack nodes.
	
	1. When a crack node is detected, an object is created and
		the constructor saves the momenta for that node.
	2.	After forces are calculated, the nodal forces are changed to map the
		initial no-contact momenta to the final contact momenta. This final
		contact momenta is calculated after updating momenta ignoring
		contact.
********************************************************************************/

#include "stdafx.h"
#include "Cracks/CrackNode.hpp"
#include "Nodes/NodalPoint.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Exceptions/CommonException.hpp"

// global point to contact conditions
vector< CrackNode * > CrackNode::crackContactNodes;

#pragma mark CrackNode: Constructors and Destructors

// Constructors
CrackNode::CrackNode(NodalPoint *nd,int flags,CrackNode *prev)
{
	theNode = nd;
	prevNode = prev;
	hasFlags = flags;
}

#pragma mark CrackNode: Methods

// check contact on this node during update strains last
void CrackNode::NodalCrackContact(double deltime,int passType)
{
	theNode->CrackContact(passType,deltime,hasFlags);
}

// next BC accessors
void CrackNode::SetPrevNode(CrackNode *next) { prevNode=next; }
CrackNode *CrackNode::GetPrevNode(void) { return prevNode; }

#pragma mark CrackNode: Class methods

// On last pass (for USAVG or SZS), will already know which
// nodes are crack nodes and now need to adjust forces
void CrackNode::ContactOnKnownNodes(double deltime,int passType)
{
	// anything to do?
	long numCrackNodes = crackContactNodes.size();
	if(numCrackNodes==0) return;
	
	// prepare for parallel
	CommonException *ccErr = NULL;
	int numNodesPerProc = (int)((double)numCrackNodes/(double)(fmobj->GetNumberOfProcessors()));
	
#pragma omp parallel if(numNodesPerProc>1)
	{
#pragma omp for
		for(long i=0;i<numCrackNodes;i++)
		{	try
			{	crackContactNodes[i]->NodalCrackContact(deltime,passType);
			}
			catch(...)
			{	if(ccErr==NULL)
				{
#pragma omp critical (error)
					ccErr = new CommonException("Unexpected error","MaterialContactNode::ContactOnKnownNodes");
				}
			}
			
		}
	}
	
	// throw error now
	if(ccErr!=NULL) throw *ccErr;
}

// delete the contact node (which delink to node), clear vector of pointers
void CrackNode::ReleaseContactNodes(void)
{
	long numCrackNodes = crackContactNodes.size();
	if(numCrackNodes==0) return;
	
	for(long i=0;i<numCrackNodes;i++)
		delete crackContactNodes[i];
	
	crackContactNodes.clear();
}

