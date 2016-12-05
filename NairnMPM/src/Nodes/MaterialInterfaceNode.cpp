/********************************************************************************
    MaterialInterfaceNode.cpp
    nairn-mpm-fea
 
    Created by John Nairn on 5/3/12.
    Copyright (c) 2102 John A. Nairn, All rights reserved.
 
    This object handles multimaterial nodes with interface laws.
 
    1. When a material node is detected, an object is created
********************************************************************************/

#include "stdafx.h"
#include "Nodes/MaterialInterfaceNode.hpp"

// global point to contact conditions
MaterialInterfaceNode *MaterialInterfaceNode::currentIntNode=NULL;

#pragma mark MaterialInterfaceNode: Constructors and Destructors

// Constructors
MaterialInterfaceNode::MaterialInterfaceNode(NodalPoint *nd,int vf,int i,int j,
						Vector *fImp,double iEnergy,MaterialInterfaceNode *last)
{
	theNode=nd;                 // the node with an interface
    vfld = vf;                  // the crack velocity field
    mati = i;                   // the material velocity field
    matipaired = j;             // the other material (or -1 if more than one other material)
	prevBC = last;
    CopyVector(&traction, fImp);
    energy = iEnergy;
}

#pragma mark MaterialInterfaceNode: Methods

// Add internal force to all nodes with imperfect interfaces 
MaterialInterfaceNode *MaterialInterfaceNode::InterfaceForceBetweenMaterials(void)
{
	theNode->MaterialInterfaceForce(this);
	return prevBC;
}

// next BC accessors
void MaterialInterfaceNode::SetPrevBC(MaterialInterfaceNode *next) { prevBC=next; }
MaterialInterfaceNode *MaterialInterfaceNode::GetPrevBC(void) { return prevBC; }

// retrieve interface calculations
double MaterialInterfaceNode::GetInterfaceTraction(Vector *fImp)
{
    CopyVector(fImp, &traction);
    return energy;
}

// retrieve velocity fields
void MaterialInterfaceNode::GetFieldInfo(int *cvf,int *i,int *j)
{   *cvf = vfld;
    *i = mati;
    *j = matipaired;
}

#pragma mark CrackNode: Class methods

// Delete all dynamically created contact BCs and restore
// currentIntNode to NULL - called in Task 0 or initialization
void MaterialInterfaceNode::RemoveInterfaceNodes(void)
{
	MaterialInterfaceNode *prevBC;
	
	while(currentIntNode!=NULL)
	{	prevBC = currentIntNode->GetPrevBC();
		delete currentIntNode;
		currentIntNode = prevBC;
	}
}

// When there are imperfect interfaces between materials, add to nodal internal force
// and track interface energy
void MaterialInterfaceNode::MaterialInterfaceOnKnownNodes(void)
{
	MaterialInterfaceNode *prevBC = currentIntNode;
	while(prevBC!=NULL)
		prevBC = prevBC->InterfaceForceBetweenMaterials();
}




