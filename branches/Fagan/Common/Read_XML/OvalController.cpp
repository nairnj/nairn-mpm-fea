/********************************************************************************
    OvalController.cpp
    NairnFEA and NairnMPM

    Created by John Nairn on 4/23/12.
    Copyright (c) 2012 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_XML/OvalController.hpp"

#pragma mark OvalController: constructors, destructor, initializers

OvalController::OvalController(int block) : ShapeController(block)
{
}

// called after initialization is done
bool OvalController::FinishSetup(void)
{   
	ShapeController::FinishSetup();
	
	// get center of the oval
	x0=(xmin+xmax)/2.;
	y0=(ymin+ymax)/2.;
	
	return TRUE;
}

#pragma mark OvalController: methods

// return true if point is in this oval centered on (x0,y0)
bool OvalController::ContainsPoint(Vector& pt)
{	double a=xmax-x0;
	double b=ymax-y0;
	return (((pt.x-x0)*(pt.x-x0)/a/a+(pt.y-y0)*(pt.y-y0)/b/b) <= 1.) ? TRUE : FALSE;
}

#pragma mark OvalController: accessors

// type of object
const char *OvalController::GetShapeName(void) { return "Oval"; }
