/********************************************************************************
    BodyPolygonController.cpp
    NairnFEA
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_MPM/BodyPolygonController.hpp"

/********************************************************************************
	BodyPolygonController: contructors and destructors
********************************************************************************/

BodyPolygonController::~BodyPolygonController()
{	xpt.clear();
	ypt.clear();
}

/********************************************************************************
	BodyPolygonController: methods
********************************************************************************/

// return true if point is in this body
bool BodyPolygonController::ContainsPoint(Vector& pt)
{
	unsigned i,crossings=0;
	double x1,x2,y1,y2,d;
	
	x1=xpt[0];
	y1=ypt[0];
	for(i=1;i<xpt.size();i++)
	{	x2=xpt[i];
		y2=ypt[i];
		d=(pt.y-y1)*(x2-x1)-(pt.x-x1)*(y2-y1);
		if((y1>=pt.y) != (y2>=pt.y))
			crossings+= (y2-y1>=0.) ? d>=0. : d<=0. ;
		if(!d && fmin(x1,x2)<=pt.x && pt.x<=fmax(x1,x2) &&
							fmin(y1,y2)<=pt.y && pt.y<=fmax(y1,y2))
			return TRUE;
		x1=x2;
		y1=y2;
	}
	return (crossings & 0x01);
}

// called after initialization is done, need to wait for pt parameters
bool BodyPolygonController::FinishSetup(void) {	return FALSE; }

// called after initialization is done
void BodyPolygonController::FinishParameter(void)
{	xpt.push_back(xparm);
	ypt.push_back(yparm);
}

/********************************************************************************
	BodyPolygonController: accessors
********************************************************************************/

// set a property
void BodyPolygonController::SetParameter(char *aName,char *value)
{
	if(strcmp(aName,"x")==0)
		sscanf(value,"%lf",&xparm);
	else if(strcmp(aName,"y")==0)
		sscanf(value,"%lf",&yparm);
}

// return if has enough parameters to use. If yes, close the  polygon
// Warning - only call once, just before use and exception if fails
bool BodyPolygonController::HasAllParameters(void)
{	if(xpt.size()<3) return FALSE;
	xpt.push_back(xpt[0]);
	ypt.push_back(ypt[0]);
	return TRUE;
}
