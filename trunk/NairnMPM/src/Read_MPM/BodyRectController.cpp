/********************************************************************************
    BodyRectController.cpp
    NairnFEA
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_MPM/BodyRectController.hpp"

/********************************************************************************
	BodyRectController: methods
********************************************************************************/

// return true if point is in this body
bool BodyRectController::ContainsPoint(Vector& pt)
{	if(pt.y>ymax || pt.y<ymin) return FALSE;
	if(pt.x>xmax || pt.x<xmin) return FALSE;
	return TRUE;
}

// type of object
const char *BodyRectController::GetShapeName(void) { return "Rectangle"; }
