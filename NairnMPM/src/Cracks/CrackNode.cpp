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

// global point to contact conditions
CrackNode *CrackNode::currentCNode = NULL;

#pragma mark CrackNode: Constructors and Destructors

// Constructors
CrackNode::CrackNode(NodalPoint *nd,CrackNode *prev)
{
	theNode = nd;
	prevNode = prev;
}

#pragma mark CrackNode: Methods

// check contact on this node during update momentum taks
CrackNode *CrackNode::NodalCrackContactAndForces(double deltime)
{
	theNode->CrackContact(UPDATE_MOMENTUM_CALL,deltime,NULL,NULL);
	return prevNode;
}

// check contact on this node during update strains last
CrackNode *CrackNode::NodalCrackContact(void)
{
	theNode->CrackContact(UPDATE_STRAINS_LAST_CALL,0.,NULL,NULL);
	return prevNode;
}

// next BC accessors
void CrackNode::SetPrevNode(CrackNode *next) { prevNode=next; }
CrackNode *CrackNode::GetPrevNode(void) { return prevNode; }

#pragma mark CrackNode: Class methods

// Delete all dynamically created contact BCs and restore
// currentCNode to NULL - called in Task 0 or initialization
void CrackNode::RemoveCrackNodes(void)
{
	CrackNode *prevCNode;
	
	while(currentCNode!=NULL)
	{	prevCNode = currentCNode->GetPrevNode();
		delete currentCNode;
		currentCNode = prevCNode;
	}
}

// In task 4, have to check if the momentum update caused new contact
// If yes, change momentum again and change total force
void CrackNode::CrackContactTask4(double deltime)
{
	CrackNode *prevCNode = currentCNode;
	while(prevCNode!=NULL)
		prevCNode = prevCNode->NodalCrackContactAndForces(deltime);
}

// On last pass (for USAVG or SZS), will already know which
// nodes are crack nodes and now need to adjust forces
void CrackNode::ContactOnKnownNodes(void)
{
	CrackNode *prevCNode = currentCNode;
	while(prevCNode!=NULL)
		prevCNode = prevCNode->NodalCrackContact();
}
