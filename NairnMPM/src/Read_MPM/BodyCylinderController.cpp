/********************************************************************************
    BodyCylinderController.hpp
    NairnFEA
    
    Created by John Nairn on 9/1/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		BodyObjectController.hpp
		BodySphereController.hpp
********************************************************************************/

#include "BodyCylinderController.hpp"

/********************************************************************************
	BodyCylinderController: constructors and destructors
********************************************************************************/

BodyCylinderController::BodyCylinderController() : BodySphereController()
{
	axis=0;
}

/********************************************************************************
	BodySphereController: methods
********************************************************************************/

// called after initialization is done
bool BodyCylinderController::FinishSetup(void)
{
	BodySphereController::FinishSetup();
	if(axis<1 || axis>3)
        ThrowSAXException("cylinder axis not set or not 1, 2, or 3");
	return TRUE;
}

// return true if point is in this body
bool BodyCylinderController::ContainsPoint(Vector& pt)
{	if(pt.y>ymax || pt.y<ymin) return FALSE;
	if(pt.x>xmax || pt.x<xmin) return FALSE;
	if(pt.z>zmax || pt.z<zmin) return FALSE;
	
	double a=xmax-x0;
	double b=ymax-y0;
	double c=zmax-z0;
	if(axis==1)
		return (((pt.y-y0)*(pt.y-y0)/b/b+(pt.z-z0)*(pt.z-z0)/c/c) <= 1.) ? TRUE : FALSE;
	else if(axis==2)
		return (((pt.x-x0)*(pt.x-x0)/a/a+(pt.z-z0)*(pt.z-z0)/c/c) <= 1.) ? TRUE : FALSE;
	return (((pt.x-x0)*(pt.x-x0)/a/a+(pt.y-y0)*(pt.y-y0)/b/b) <= 1.) ? TRUE : FALSE;
}

/********************************************************************************
	BodyCylinderController: accessors
********************************************************************************/

// set a property
void BodyCylinderController::SetBodyProperty(const char *aName,char *value,CommonReadHandler *reader)
{
	if(strcmp(aName,"axis")==0)
	{	sscanf(value,"%d",&axis);
	}
	else
		BodyObjectController::SetBodyProperty(aName,value,reader);
}

// type of object
const char *BodyCylinderController::GetObjectType(void) { return "Cylinder"; }

