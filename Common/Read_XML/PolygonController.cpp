/********************************************************************************
    PolygonController.cpp
    nairn-mpm-fea
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_XML/PolygonController.hpp"

#pragma mark PolygonController: initializers

// constructor
PolygonController::PolygonController(int block) : ShapeController(block)
{
}

// destructor
PolygonController::~PolygonController()
{	xpt.clear();
	ypt.clear();
    xparm = 0.;
    yparm = 0.;
}

// set (x,y) point
void PolygonController::SetParameter(const char *aName,const char *value)
{
	if(strcmp(aName,"x")==0)
	{	sscanf(value,"%lf",&xparm);
		xparm*=distScaling;
	}
	else if(strcmp(aName,"y")==0)
	{	sscanf(value,"%lf",&yparm);
		yparm*=distScaling;
	}
}

// called after initialization is done
// not thread safe due to push_back()
bool PolygonController::FinishParameter(void)
{	xpt.push_back(xparm);
	ypt.push_back(yparm);
	return true;
}

// add point using numbers
void PolygonController::AddPoint(double xp,double yp)
{	xpt.push_back(xp);
	ypt.push_back(yp);
}

// called after initialization is done, need to wait for pt parameters
bool PolygonController::FinishSetup(void) {	return FALSE; }

// return if has enough parameters to use. If yes, close the  polygon
// Warning - only call once, just before use and exception if fails
// not thread safe due to push_back()
bool PolygonController::HasAllParameters(void)
{	if(xpt.size()<3) return FALSE;
	xpt.push_back(xpt[0]);
	ypt.push_back(ypt[0]);
	return TRUE;
}

#pragma mark PolygonController: methods

// return true if point is in this body
bool PolygonController::ContainsPoint(Vector &pt)
{
	unsigned i,crossings=0;
	double x1,x2,y1,y2,d;
	
	x1=xpt[0];
	y1=ypt[0];
	for(i=1;i<xpt.size();i++)
	{	x2=xpt[i];
		y2=ypt[i];
		d=(pt.y-y1)*(x2-x1)-(pt.x-x1)*(y2-y1);
		
		// get crossing unless both y's on same side of edge
		if((y1>=pt.y) != (y2>=pt.y))
			crossings+= (y2-y1>=0.) ? d>=0. : d<=0. ;
		
		// if d is 0, check if point is on line (and thus in polygon)
		if(!d && fmin(x1,x2)<=pt.x && pt.x<=fmax(x1,x2) &&
							fmin(y1,y2)<=pt.y && pt.y<=fmax(y1,y2))
			return TRUE;
		x1=x2;
		y1=y2;
	}
	return (crossings & 0x01);
}

#pragma mark PolygonController: accessors

// type of object
const char *PolygonController::GetShapeName(void) { return "Polygon"; }





