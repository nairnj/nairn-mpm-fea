/********************************************************************************
    BodyBoxController.cpp
    NairnFEA
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_MPM/BodyBoxController.hpp"

/********************************************************************************
	BodyBoxController: methods
********************************************************************************/

// return true if point is in this body
bool BodyBoxController::ContainsPoint(Vector& pt)
{	if(pt.y>ymax || pt.y<ymin) return FALSE;
	if(pt.x>xmax || pt.x<xmin) return FALSE;
	if(pt.z>zmax || pt.z<zmin) return FALSE;
	return TRUE;
}

/********************************************************************************
	BodyBoxController: accessors
********************************************************************************/

// override for 3D objects
bool BodyBoxController::Is2DBodyObject(void) { return FALSE; }

// type of object
char *BodyBoxController::GetObjectType(void) { return "Box"; }

