/********************************************************************************
    ArcController.cpp
    NairnFEA
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_XML/ArcController.hpp"

/********************************************************************************
	ArcController: Constructors and Destructor
********************************************************************************/

ArcController::ArcController(int block) : LineController(block,true)
{	startAngle=0.;
	endAngle=2.*PI_CONSTANT;
}

// called after initialization is done to precalculate terms helpful in ContainsPoint
bool ArcController::FinishSetup(void)
{
    ShapeController::FinishSetup();
	centerX=(xmax+xmin)/2.0;
	centerY=(ymax+ymin)/2.0;
	a=(xmax-xmin)/2.0;
	b=(ymax-ymin)/2.0;
    return true;
}

/********************************************************************************
	ArcController: methods
********************************************************************************/

// Determine if point is sufficiently close to the line from
//   (xmin,ymin) to (xmax,ymax). Must be within rectangle a
//   distance tolerance from the line in all directions 
bool ArcController::ContainsPoint(Vector& v)
{
	double deltaX, deltaY, sita=0.0, R;

	deltaX=v.x-centerX;
	deltaY=v.y-centerY;

	// convert to polar coordinates
	double polarR=sqrt(deltaX*deltaX+deltaY*deltaY);
	if(polarR>0)
	{	sita=acos(deltaX/polarR);			// 0 to PI
		if(deltaY<0) sita=2.0*PI_CONSTANT-sita;
	}
	
	if(a==b)
		R=a;
	else
		R=a*b/sqrt(a*a*sin(sita)*sin(sita)+b*b*cos(sita)*cos(sita));
	
	// get positive angle only (0 to 2PI)
	if((sita>=startAngle) && (sita<=endAngle) && (fabs(polarR-R)<=tolerance)) return true;
	
	// does arc have negative angles
	if(startAngle>=0. || deltaY>=0.) return false;
	
	// check the negative angle
	sita -= 2.0*PI_CONSTANT;
	if((sita>=startAngle) && (sita<=endAngle) && (fabs(polarR-R)<=tolerance)) return true;

    return false;
}

// set a property
void ArcController::SetProperty(const char *aName,char *value,CommonReadHandler *reader)
{
	if(strcmp(aName,"start")==0)
	{	sscanf(value,"%lf",&startAngle);
		startAngle*=PI_CONSTANT/180.;		// convert to radians
	}
	else if(strcmp(aName,"end")==0)
	{	sscanf(value,"%lf",&endAngle);
		endAngle*=PI_CONSTANT/180.;		// convert to radians
	}
	else
		LineController::SetProperty(aName,value,reader);
}

// type of object
const char *ArcController::GetShapeName(void) { return "Arc"; }


