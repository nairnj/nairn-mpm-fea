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
	BodyObjectController: constructors and destructors
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

/********************************************************************************
	BodyObjectController: methods
********************************************************************************/

// return true if point is in this body
bool BodyObjectController::ContainsPoint(Vector& pt) { return FALSE; }

// called after initialization is done, return TRUE if ready to use
// or FALSE if this object needs to wait for parameters
bool BodyObjectController::FinishSetup(void)
{
    double temp;
    if(xmin>xmax)
	{	temp=xmax;
        xmax=xmin;
        xmin=temp;
    }
    if(DbleEqual(xmin,xmax))
        ThrowSAXException("%s: xmax cannot equal xmin in input parameters.",GetObjectType());
		
    if(ymin>ymax)
	{	temp=ymax;
        ymax=ymin;
        ymin=temp;
    }
    if(DbleEqual(ymin,ymax))
        ThrowSAXException("%s: ymax cannot equal ymin in input parameters.",GetObjectType());
	
	if(!Is2DBodyObject())
	{	if(zmin>zmax)
		{	temp=zmax;
			zmax=zmin;
			zmin=temp;
		}
		if(DbleEqual(zmin,zmax))
			ThrowSAXException("%s: zmax cannot equal zmin in input parameters.",GetObjectType());
	}
	
	return TRUE;
}

// called after one parameter setting is done
bool BodyObjectController::FinishParameter(void) { return TRUE; }

/********************************************************************************
	BodyObjectController: accessors
********************************************************************************/

// set a property
void BodyObjectController::SetBodyProperty(const char *aName,char *value,CommonReadHandler *reader)
{
	if(strcmp(aName,"xmin")==0)
		xmin=reader->ReadX(value,distScaling);
	else if(strcmp(aName,"ymin")==0)
		ymin=reader->ReadY(value,distScaling);
	else if(strcmp(aName,"zmin")==0)
		zmin=reader->ReadZ(value,distScaling);
	else if(strcmp(aName,"xmax")==0)
		xmax=reader->ReadX(value,distScaling);
	else if(strcmp(aName,"ymax")==0)
		ymax=reader->ReadY(value,distScaling);
	else if(strcmp(aName,"zmax")==0)
		zmax=reader->ReadZ(value,distScaling);
}

// set an object parameter (in subordinate command)
void BodyObjectController::SetParameter(char *aName,char *value) { }

// override for 3D objects
bool BodyObjectController::Is2DBodyObject(void) { return TRUE; }

// set the scaling
void BodyObjectController::SetScaling(double scale) { distScaling=scale; }

// called at most once right before use to see if enough parameters have been set
bool BodyObjectController::HasAllParameters(void) { return TRUE; }

// to allow obbject to decode object-specific character data
bool BodyObjectController::SetBodyPropertyFromData(char *bData,CommonReadHandler *reader)
{	return TRUE;
}

// type of object
char *BodyObjectController::GetObjectType(void) { return "Body"; }
