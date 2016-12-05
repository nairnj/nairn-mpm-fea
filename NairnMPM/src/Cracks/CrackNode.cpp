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
	prevBC = prev;
}

#pragma mark CrackNode: Methods

// Add internal force to all nodes with imperfect interfaces 
CrackNode *CrackNode::InterfaceForceOnCrack(void)
{
	theNode->CrackInterfaceForce();
	return prevBC;
}

// check contact on this node during update momentum taks
CrackNode *CrackNode::NodalCrackContactAndForces(double deltime)
{
	theNode->CrackContact(TRUE,deltime,NULL,NULL);
	return prevBC;
}

// check contact on this node during update strains last
CrackNode *CrackNode::NodalCrackContact(void)
{
	theNode->CrackContact(FALSE,0.,NULL,NULL);
	return prevBC;
}

// next BC accessors
void CrackNode::SetPrevBC(CrackNode *next) { prevBC=next; }
CrackNode *CrackNode::GetPrevBC(void) { return prevBC; }

#pragma mark CrackNode: Class methods

// Delete all dynamically created contact BCs and restore
// currentCNode to NULL - called in Task 0 or initialization
void CrackNode::RemoveCrackNodes(void)
{
	CrackNode *prevBC;
	
	while(currentCNode!=NULL)
	{	prevBC = currentCNode->GetPrevBC();
		delete currentCNode;
		currentCNode = prevBC;
	}
}

// In task 4, have to check if the momentum update caused new contact
// If yes, change momentum again and change total force
void CrackNode::CrackContactTask4(double deltime)
{
	CrackNode *prevBC = currentCNode;
	while(prevBC!=NULL)
		prevBC = prevBC->NodalCrackContactAndForces(deltime);
}

// On last pass (for USAVG or SZS), will already know which
// nodes are crack nodes and now need to adjust forces
void CrackNode::ContactOnKnownNodes(void)
{
	CrackNode *prevBC = currentCNode;
	while(prevBC!=NULL)
		prevBC = prevBC->NodalCrackContact();
}
	
// When there are imperfect interfaces, add to nodal internal force
void CrackNode::CrackInterfaceOnKnownNodes(void)
{
	CrackNode *prevBC = currentCNode;
	while(prevBC!=NULL)
		prevBC = prevBC->InterfaceForceOnCrack();
}

	


