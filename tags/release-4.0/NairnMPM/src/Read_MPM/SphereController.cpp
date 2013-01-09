/********************************************************************************
    SphereController.cpp
    NairnMPM
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_MPM/SphereController.hpp"

#pragma mark SphereController: initializers

// constructor
SphereController::SphereController(int block) : ShapeController(block)
{
}

// called after initialization is done
bool SphereController::FinishSetup(void)
{
	ShapeController::FinishSetup();
	
	// get center of the oval
	x0=(xmin+xmax)/2.;
	y0=(ymin+ymax)/2.;
	z0=(zmin+zmax)/2.0;
	
	return TRUE;
}

#pragma mark SphereController: methods

// return true if point is in this sphere
bool SphereController::ContainsPoint(Vector& pt)
{	double a=xmax-x0;
	double b=ymax-y0;
	double c=zmax-z0;
	return (((pt.x-x0)*(pt.x-x0)/a/a+(pt.y-y0)*(pt.y-y0)/b/b+(pt.z-z0)*(pt.z-z0)/c/c) <= 1.) ? TRUE : FALSE;
}

#pragma mark SphereController: accessors

// override for 3D objects
bool SphereController::Is2DShape(void) { return FALSE; }

// type of object
const char *SphereController::GetShapeName(void) { return "Sphere"; }

