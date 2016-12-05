/********************************************************************************
    RectController.cpp
    NairnFEA
    
    Created by John Nairn on 8/8/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_XML/RectController.hpp"

#pragma mark RectController: constructors, destructor, initializers

RectController::RectController(int block) : ShapeController(block)
{
	startAngle = -1.;
}

// set (x,y) point
void RectController::SetParameter(const char *aName,const char *value)
{
	if(strcmp(aName,"start")==0)
	{	sscanf(value,"%lf",&startAngle);
	}
	else if(strcmp(aName,"end")==0)
	{	sscanf(value,"%lf",&endAngle);
	}
}

// called after initialization is done
bool RectController::FinishParameter(void)
{	if(startAngle<0. || startAngle>360.) return false;
	if(endAngle < startAngle) return false;
	
	// convert to radians
	startAngle *= PI_CONSTANT/180.;
	endAngle *= PI_CONSTANT/180.;
	
	return true;
}

#pragma mark RectController: methods

// return true if point is in this rectangle
bool RectController::ContainsPoint(Vector& v)
{	if(v.x<=xmax && v.x>=xmin && v.y<=ymax && v.y>=ymin)
	{	return (startAngle>=0.) ? CheckAngle(v.x-0.5*(xmax+xmin),v.y-0.5*(ymax+ymin)) : true;
	}
	return false;
}

// get angle from 0 to 2*pi (assume startAngle>=0)
bool RectController::CheckAngle(double x,double y)
{	double len = sqrt(x*x+y*y);
	if(DbleEqual(len,0.)) return true;
	double angle;
	if(x>0.)
		angle = (y>0.) ? asin(y/len) : asin(y/len)+2*PI_CONSTANT;
	else
		angle = PI_CONSTANT-asin(y/len);
	if(angle > endAngle) return false;
	if(angle < startAngle)
	{	if(endAngle<=2*PI_CONSTANT) return false;
		// in case endAngle > 360
		if(angle>endAngle-2*PI_CONSTANT) return false;
	}
	return true;
}

#pragma mark RectController: accessors

// type of object
const char *RectController::GetShapeName(void) { return "Rectangle"; }

