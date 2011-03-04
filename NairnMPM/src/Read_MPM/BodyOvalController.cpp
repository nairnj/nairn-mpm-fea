/********************************************************************************
    BodyOvalController.cpp
    NairnFEA
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_MPM/BodyOvalController.hpp"

/********************************************************************************
	BodyOvalController: methods
********************************************************************************/

// called after initialization is done
bool BodyOvalController::FinishSetup(void)
{
	BodyObjectController::FinishSetup();
	
	// get center of the oval
	x0=(xmin+xmax)/2.;
	y0=(ymin+ymax)/2.;
	
	return TRUE;
}

// return true if point is in this body
bool BodyOvalController::ContainsPoint(Vector& pt)
{	double a=xmax-x0;
	double b=ymax-y0;
	return (((pt.x-x0)*(pt.x-x0)/a/a+(pt.y-y0)*(pt.y-y0)/b/b) <= 1.) ? TRUE : FALSE;
}

// type of object
char *BodyOvalController::GetObjectType(void) { return "Oval"; }
