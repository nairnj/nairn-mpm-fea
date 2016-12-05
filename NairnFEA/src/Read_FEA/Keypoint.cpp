/********************************************************************************
    Keypoint.hpp
    NairnFEA
    
    Created by John Nairn on 6/23/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_FEA/Keypoint.hpp"
#include "Read_XML/NodesController.hpp"

/********************************************************************************
	Keypoins: Constructors and Destructor
********************************************************************************/

// Initialize, but assumes name already 32 characters or less and unique
// ... which is done by prior call to keyPts->ValidName()
// ... but internally created keypoints will have empty name
Keypoint::Keypoint(const char *keyName,double xpt,double ypt)
{
	strcpy(keyID,keyName);
	x=xpt;
	y=ypt;
	node=-1;
	temperature=0.0;
}

/********************************************************************************
	Keypoint: methods
********************************************************************************/

// Create node for this keypoint
void Keypoint::AssignToNode(void)
{
	if(node>=0) return;
	theNodes->AddNode(x,y,(double)0.,temperature);
	node=theNodes->numObjects;
}

// are two keypoints at the same coordinates
int Keypoint::SamePtAs(Keypoint *otherKey)
{
	Vector keyCoords;
	otherKey->GetKeypointXY(&keyCoords);
	if(!DbleEqual(x,keyCoords.x)) return FALSE;
	if(!DbleEqual(y,keyCoords.y)) return FALSE;
	return TRUE;
}

/********************************************************************************
	Keypoint: accessors
********************************************************************************/

void Keypoint::GetKeypointXY(Vector *keyCoords)
{	keyCoords->x=x;
	keyCoords->y=y;
}

