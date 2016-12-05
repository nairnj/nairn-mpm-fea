/********************************************************************************
    OvalController.cpp
    nairn-mpm-fea

    Created by John Nairn on 4/23/12.
    Copyright (c) 2012 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_XML/OvalController.hpp"

#pragma mark OvalController: constructors, destructor, initializers

OvalController::OvalController(int block) : RectController(block)
{
}

// called after initialization is done
bool OvalController::FinishSetup(void)
{   
	ShapeController::FinishSetup();
	
	// get center of the oval
	x0=(xmin+xmax)/2.;
	y0=(ymin+ymax)/2.;
	
	// x and y direction radii
	a=xmax-x0;
	b=ymax-y0;
	
	return TRUE;
}

#pragma mark OvalController: methods

// return true if point is in this oval centered on (x0,y0)
bool OvalController::ContainsPoint(Vector& pt)
{	double x = pt.x-x0;
	double y = pt.y-y0;
	
	// exit if not in ellipse: (x/a)^2+(y/b)^2 = 1
	if((x*x/(a*a)+y*y/(b*b)) > 1.) return false;
	
	// within an arc?
	return startAngle>=0. ? CheckAngle(x,y) : true;
}

#pragma mark OvalController: accessors

// type of object
const char *OvalController::GetShapeName(void) { return "Oval"; }
