/********************************************************************************
    BodySphereController.cpp
    NairnFEA
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_MPM/BodySphereController.hpp"

/********************************************************************************
	BodySphereController: methods
********************************************************************************/

// called after initialization is done
bool BodySphereController::FinishSetup(void)
{
	BodyObjectController::FinishSetup();
	
	// get center of the oval
	x0=(xmin+xmax)/2.;
	y0=(ymin+ymax)/2.;
	z0=(zmin+zmax)/2.0;
	
	return TRUE;
}

// return true if point is in this body
bool BodySphereController::ContainsPoint(Vector& pt)
{	double a=xmax-x0;
	double b=ymax-y0;
	double c=zmax-z0;
	return (((pt.x-x0)*(pt.x-x0)/a/a+(pt.y-y0)*(pt.y-y0)/b/b+(pt.z-z0)*(pt.z-z0)/c/c) <= 1.) ? TRUE : FALSE;
}

/********************************************************************************
	BodyBoxController: accessors
********************************************************************************/

// override for 3D objects
bool BodySphereController::Is2DShape(void) { return FALSE; }

// type of object
const char *BodySphereController::GetShapeName(void) { return "Sphere"; }

