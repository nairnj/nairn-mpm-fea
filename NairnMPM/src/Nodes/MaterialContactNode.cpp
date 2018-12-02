/********************************************************************************
	MaterialContactNode.cpp
	nairn-mpm-fea
	 
	Created by John Nairn on 3/15/2018.
	Copyright (c) 2018 John A. Nairn, All rights reserved.
		
	This object handles material contact nodes when using machine learning
	method to find the normals
********************************************************************************/

#include "stdafx.h"
#include "Nodes/MaterialContactNode.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Exceptions/CommonException.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Patches/GridPatch.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Cracks/CrackSurfaceContact.hpp"

#define PARALLEL_LINKING

// global point to contact conditions
vector< MaterialContactNode * > MaterialContactNode::materialContactNodes;

#pragma mark MaterialContactNode: Constructors and Destructors

// Constructors
MaterialContactNode::MaterialContactNode(NodalPoint *nd,MaterialContactNode *prev)
{
	theNode = nd;
	prevNode = prev;
	lists = NULL;
}

// Destructor
MaterialContactNode::~MaterialContactNode()
{
//#define DEBUG_LISTS
#ifdef DEBUG_LISTS
	if(lists!=NULL)
	{	for(int i=0;i<maxCrackFields;i++)
		{	cout << "# " <<fmobj->mstep << ":lists[" << i << "](size:" << lists[i].size() << ",node:" << theNode->num << ") ";
			for(int j=0;j<lists[i].size();j++)
			{	cout << lists[i][j] << " ";
			}
			cout << endl;
		}
	}
#endif
	
	// delete the lists
	if(lists!=NULL) delete [] lists;
	
	// null pointer from node to this contact node
	theNode->contactData = NULL;
}

// Get the node
NodalPoint *MaterialContactNode::GetTheNode(void) { return theNode; }

#pragma mark MaterialContactNode: Methods

// check contact on this node during update strains last
void MaterialContactNode::NodalMaterialContact(double dtime,int passType)
{
	theNode->MaterialContactOnNode(dtime,passType,this);
}

// Add material point to list of materials points in a crack velocity field
void MaterialContactNode::AddMaterialContactPoint(int pnum,short vfld)
{
#ifdef PARALLEL_LINKING
	if(theNode->NeedsParticleListed(vfld))
	{
#pragma omp critical (pushlist)
		{
			lists[vfld].push_back(pnum);
		}
	}
#else
	if(theNode->NeedsParticleListed(vfld))
		lists[vfld].push_back(pnum);
#endif
	
}

// create lists needed to store material points on target not
void MaterialContactNode::PrepareForLists(void)
{
	// need maxCrackFields list (at most)
	//    (1 if no cracks or MAX_FIELDS_FOR_CRACKS/MAX_FIELDS_FOR_ONE_CRACK if any cracks)
	lists = new vector< int >[maxCrackFields];
}

#pragma mark MaterialContactNode: Accessors

// next node accessors
void MaterialContactNode::SetPrevNode(MaterialContactNode *next) { prevNode=next; }
MaterialContactNode *MaterialContactNode::GetPrevNode(void) { return prevNode; }

// list of particle in one velocity field
vector<int> MaterialContactNode::ParticleLists(int vfld) { return lists[vfld]; }

#pragma mark CrackNode: Class methods

// On last pass (for USAVG or SZS), will already know which
// nodes are crack nodes and now need to adjust forces
void MaterialContactNode::ContactOnKnownNodes(double dtime,int passType)
{
	// anythiing to do?
	long numContactNodes = materialContactNodes.size();
	if(numContactNodes==0) return;
	
	// prepare for parallel
	CommonException *mcErr = NULL;
	int numNodesPerProc = (int)((double)numContactNodes/(double)(fmobj->GetNumberOfProcessors()));
	
	// Do contact on all material contact nodes
#pragma omp parallel if(numNodesPerProc>1)
	{
#pragma omp for
		for(long i=0;i<numContactNodes;i++)
		{	try
			{	materialContactNodes[i]->NodalMaterialContact(dtime,passType);
			}
			catch(...)
			{	if(mcErr==NULL)
				{
#pragma omp critical (error)
					mcErr = new CommonException("Unexpected error","MaterialContactNode::ContactOnKnownNodes");
				}
			}
		}
	}
	
	// throw error now
	if(mcErr!=NULL) throw *mcErr;
}

// delete the contact node (which delink to node), clear vector of pointers
void MaterialContactNode::ReleaseContactNodes(void)
{
	long numContactNodes = materialContactNodes.size();
	if(numContactNodes==0) return;
	
	for(long i=0;i<numContactNodes;i++)
		delete materialContactNodes[i];
	
	materialContactNodes.clear();

}
