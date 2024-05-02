/********************************************************************************
    NodesController.cpp
    NairnFEA
    
    Created by John Nairn on 6/22/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Read_XML/NodesController.hpp"
#include "Nodes/NodalPoint.hpp"
#ifdef MPM_CODE
	#include "NairnMPM_Class/NairnMPM.hpp"
#else
	#include "NairnFEA_Class/NairnFEA.hpp"
	#include <algorithm>
#endif

NodesController *theNodes=NULL;

/********************************************************************************
	NodesController: methods
********************************************************************************/

#ifdef MPM_CODE

// MPM does not need temperature
void NodesController::AddNode(double x,double y,double z)
{
	NodalPoint *newNode;
	if(fmobj->IsThreeD())
		newNode = NodalPoint::Create3DNode(numObjects+1,x,y,z);
	else
		newNode = NodalPoint::Create2DNode(numObjects+1,x,y);
	AddObject(newNode);
}

#else

// add new keypoint - assumes keyName is valid and unique
// throws std::bad_alloc
void NodesController::AddNode(double x,double y,double z,double temp)
{
	NodalPoint *newNode;
	if(fmobj->IsThreeD())
		newNode = NodalPoint::Create3DNode(numObjects+1,x,y,z);
	else
		newNode = NodalPoint::Create2DNode(numObjects+1,x,y);
	AddObject(newNode);
	newNode->gTemperature = temp;
}

#endif

// assemble into array used in the code
NodalPoint **NodesController::SetNodeArray(int &countNodes,double *xmin,double *xmax,double *ymin,double *ymax,double *zmin,double *zmax)
{
	// make 1-based array of nodal points
	if(numObjects==0) return NULL;
	NodalPoint **ndArray = new (std::nothrow) NodalPoint *[numObjects+1];
	if(ndArray==NULL) return NULL;
	
	// fill the array
	NodalPoint *aNode=(NodalPoint *)firstObject;
	if(aNode!=NULL)
	{	*xmin=*xmax=aNode->x;
		*ymin=*ymax=aNode->y;
		*zmin=*zmax=aNode->z;
	}
	countNodes=0;
	while(aNode!=NULL)
	{	countNodes++;
		ndArray[countNodes]=aNode;
		if(aNode->x<*xmin) *xmin=aNode->x;
		if(aNode->x>*xmax) *xmax=aNode->x;
		if(aNode->y<*ymin) *ymin=aNode->y;
		if(aNode->y>*ymax) *ymax=aNode->y;
		if(aNode->z<*zmin) *zmin=aNode->z;
		if(aNode->z>*zmax) *zmax=aNode->z;
		aNode=(NodalPoint *)aNode->GetNextObject();
	}
	return ndArray;
}

// assemble into array used in the code when they are resequenced
NodalPoint **NodesController::SetNodeArray(int *revMap)
{
	// make 1-base array of nodal points
	if(numObjects==0) return NULL;
	NodalPoint **ndArray = new (std::nothrow) NodalPoint *[numObjects+1];
	if(ndArray==NULL) return NULL;
	
	// fill the array
	NodalPoint *aNode=(NodalPoint *)firstObject;
	nnodes=0;
	while(aNode!=NULL)
	{	nnodes++;
		ndArray[revMap[nnodes]]=aNode;
		aNode->num=revMap[nnodes];
		aNode=(NodalPoint *)aNode->GetNextObject();
	}
	return ndArray;
}

// next node number
int NodesController::NextNodeNumber(void) { return numObjects+1; }

#ifdef FEA_CODE
void NodesController::MidPoint(int *eNode,int numNodes,Vector *midPt)
{
	int sNode[MaxElNd];
	int i;
	
	// sort copy of the nodes
	for(i=0;i<numNodes;i++)
		sNode[i]=eNode[i];
	std::sort(sNode,sNode+numNodes);
	
	// scan through the nodes
	int iNode=0;
	NodalPoint *aNode=(NodalPoint *)firstObject;
	midPt->x=0.;
	midPt->y=0.;
	while(iNode<numNodes)
	{	if(aNode->num==sNode[iNode])
		{	midPt->x+=aNode->x;
			midPt->y+=aNode->y;
			iNode++;
		}
		aNode=(NodalPoint *)aNode->GetNextObject();
		if(aNode==NULL) break;
	}
	
	// get mean
	if(iNode>0)
	{	midPt->x/=(double)iNode;
		midPt->y/=(double)iNode;
	}
}
#endif
