/********************************************************************************
	ShellController.cpp
	nairn-mpm-fea

	Created by John Nairn on 1/22/2014.
	Copyright (c) 2014 John A. Nairn, All rights reserved.
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Read_MPM/ShellController.hpp"
#include "Read_XML/Expression.hpp"

#pragma mark ShellController: initializers

// constructor
ShellController::ShellController(int block) : BoxController(block)
{
}

// set a property
// throws std::bad_alloc, SAXException()
void ShellController::SetProperty(const char *aName,char *value,CommonReadHandler *reader)
{
    if(strcmp(aName,"radius")==0)
	{	// bad if NULL
		if(aName==NULL)
			ThrowSAXException("Shell radius function is missing");
		radiusFunction = Expression::CreateExpression(value,"Shell radius function of height (h) time is not valid");
	}
    else if(strcmp(aName,"thickness")==0)
	{	// bad if NULL
		if(aName==NULL)
			ThrowSAXException("Shell thickness function is missing");
		thicknessFunction = Expression::CreateExpression(value,"Shell thickness function of height (h) time is not valid");
	}
    else
	{	// pass on to get axis and axis ranges
        BoxController::SetProperty(aName,value,reader);
	}
}

// called after initialization is done
bool ShellController::FinishSetup(void) { return TRUE; }

#pragma mark SphereController: methods

// return true if point is in this sphere
bool ShellController::ContainsPoint(Vector& v)
{
	double varHeight;
	double dx,dy;
	
    if(axis==1)
    {   if(v.x>xmax || v.x<xmin) return FALSE;
		varHeight = v.x/distScaling;
		dx = v.y - ymin;
		dy = v.z - zmin;
    }
    else if(axis==2)
    {   if(v.y>ymax || v.y<ymin) return FALSE;
		varHeight = v.y/distScaling;
		dx = v.x - xmin;
		dy = v.z - zmin;
    }
    else
    {   if(v.z>zmax || v.z<zmin) return FALSE;
		varHeight = v.z/distScaling;
		dx = v.x - xmin;
		dy = v.y - ymin;
    }
	
	// find radius
	double r = sqrt(dx*dx+dy*dy);

	unordered_map<string,double> vars;
	vars["h"] = varHeight;
	
	// inner radius
	double rmin = distScaling*radiusFunction->EvaluateFunction(vars);
	if(r<rmin) return false;
	
	// thickness
	double rmax = rmin + distScaling*thicknessFunction->EvaluateFunction(vars);
	if(r>rmax) return false;
	
	// OK
	return true;
}

#pragma mark SphereController: accessors

// override for 3D objects
bool ShellController::Is2DShape(void) { return FALSE; }

// type of object
const char *ShellController::GetShapeName(void) { return "Shell"; }

