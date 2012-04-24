/********************************************************************************
    BodyObjectController.cpp
    NairnFEA
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_MPM/BodyObjectController.hpp"
#include "Read_XML/CommonReadHandler.hpp"

BodyObjectController *theBody=NULL;

/********************************************************************************
	BodyObjectController: constructors and destructors and initializers
********************************************************************************/

BodyObjectController::BodyObjectController()
{
	xmin=0.;
	xmax=0.;
	ymin=0.;
	ymax=0.;
	zmin=0.;
	zmax=0.;
	distScaling=1.;
}

BodyObjectController::~BodyObjectController() { }

// set a property
void BodyObjectController::SetProperty(const char *aName,char *value,CommonReadHandler *reader)
{
    if(strcmp(aName,"x1")==0 || strcmp(aName,"xmin")==0)
    {	xmin=reader->ReadX(value,distScaling);
    }
    else if(strcmp(aName,"y1")==0 || strcmp(aName,"ymin")==0)
    {	ymin=reader->ReadY(value,distScaling);
    }
    else if(strcmp(aName,"z1")==0 || strcmp(aName,"zmin")==0)
    {	zmin=reader->ReadZ(value,distScaling);
    }
    else if(strcmp(aName,"x2")==0 || strcmp(aName,"xmax")==0)
    {	xmax=reader->ReadX(value,distScaling);
    }
    else if(strcmp(aName,"y2")==0 || strcmp(aName,"ymax")==0)
    {	ymax=reader->ReadY(value,distScaling);
    }
    else if(strcmp(aName,"z2")==0 || strcmp(aName,"zmax")==0)
    {	zmax=reader->ReadZ(value,distScaling);
    }
}

// set a property from value without read handler
void BodyObjectController::SetProperty(const char *aName,double value)
{
	if(strcmp(aName,"x1")==0 || strcmp(aName,"xmin")==0)
		xmin=value*distScaling;
	else if(strcmp(aName,"y1")==0 || strcmp(aName,"ymin")==0)
		ymin=value*distScaling;
	else if(strcmp(aName,"z1")==0 || strcmp(aName,"zmin")==0)
		zmin=value*distScaling;
	else if(strcmp(aName,"x2")==0 || strcmp(aName,"xmax")==0)
		xmax=value*distScaling;
	else if(strcmp(aName,"y2")==0 || strcmp(aName,"ymax")==0)
		ymax=value*distScaling;
	else if(strcmp(aName,"z2")==0 || strcmp(aName,"zmax")==0)
		zmax=value*distScaling;
}

// to allow object to decode object-specific character data
// throw an exception if bad data
void BodyObjectController::SetProperty(char *bData,CommonReadHandler *reader) {}

// set the scaling
void BodyObjectController::SetScaling(double scale) { distScaling=scale; }

// set an object parameter (in subordinate command)
// called for attributes on XML objects subordinate to the shape command
void BodyObjectController::SetParameter(const char *aName,const char *value) { }

// called after finish attributes of subordinate command
// return FALSE if not set correctly, or TRUE is OK to continue
bool BodyObjectController::FinishParameter(void) { return TRUE; }

// called after initialization is done, return TRUE if ready to use
// or FALSE if this object needs to wait for parameters
// This base class requires min and max (x, y and z) to differ and
//      reorders if needed. This it correct for rect, oval, box, sphere
//      and cylinder, but maybe not for others.
bool BodyObjectController::FinishSetup(void)
{
    double temp;
    if(xmin>xmax)
	{	temp=xmax;
        xmax=xmin;
        xmin=temp;
    }
    if(DbleEqual(xmin,xmax))
        ThrowSAXException("%s: xmax cannot equal xmin in input parameters.",GetShapeName());
    
    if(ymin>ymax)
	{	temp=ymax;
        ymax=ymin;
        ymin=temp;
    }
    if(DbleEqual(ymin,ymax))
        ThrowSAXException("%s: ymax cannot equal ymin in input parameters.",GetShapeName());
	
	if(!Is2DShape())
	{	if(zmin>zmax)
    {	temp=zmax;
        zmax=zmin;
        zmin=temp;
    }
		if(DbleEqual(zmin,zmax))
			ThrowSAXException("%s: zmax cannot equal zmin in input parameters.",GetShapeName());
	}
	
	return TRUE;
}

// some shapes might call this right be fore use. Return TRUE or FALSE
// if has all parameters. Normally only for shapes with subordinate commands.
bool BodyObjectController::HasAllParameters(void) { return TRUE; }

/********************************************************************************
	BodyObjectController: methods
********************************************************************************/

// return true if point is in this body
bool BodyObjectController::ContainsPoint(Vector& pt) { return FALSE; }

/********************************************************************************
	BodyObjectController: accessors
********************************************************************************/

// override for 3D objects
bool BodyObjectController::Is2DShape(void) { return TRUE; }

// type of object
const char *BodyObjectController::GetShapeName(void) { return "Shape"; }

