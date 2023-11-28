/********************************************************************************
    LineController.cpp
    NairnFEA
    
    Created by John Nairn on 6/62/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Elements/ElementBase.hpp"
#include "Read_XML/LineController.hpp"

#pragma mark LineController: Constructors and Destructor

// Create, may be 2D or 3D line
LineController::LineController(int block,bool is2D) : ShapeController(block)
{
	tolerance = 0.;
	twoDShape = is2D;
}

// used in FEA and always 2D
LineController::LineController(int block,double x1,double x2,double y1,double y2,double tolerate)
		: ShapeController(block,x1,x2,y1,y2)
{
	// globals for length and tolerance
	tolerance = tolerate;
	distanceSq=(xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin);
	dtolerance=sqrt(distanceSq)*tolerate;
}

// called after initialization is done
bool LineController::FinishSetup(void)
{
	// globals for length and tolerance
	if(twoDShape)
		distanceSq=(xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin);
	else
		distanceSq=(xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)+(zmax-zmin)*(zmax-zmin);
	dtolerance=sqrt(distanceSq)*tolerance;
    return TRUE;
}

#pragma mark LineController: methods

// Deterime if point is sufficiently close to the line from
//   (xmin,ymin) to (xmax,ymax). Must be within rectangle a
//   distance tolerance from the line in all directions 
bool LineController::ContainsPoint(Vector& v)
{
	double dx = xmax-xmin, dy = ymax-ymin;
	double px = v.x-xmin, py = v.y-ymin;
	
	if(twoDShape)
	{	// Let t = (-dy,dx)/len be tangent unit vector, then (t.p)*len
		// is distance tangent to line*len. Compare to tolerance
		double ypr = -dy*px + dx*py;
 		if(ypr>dtolerance) return false;
		if(ypr<-dtolerance) return false;
		
		// Let n = (dx,xy)/len be unit vector along the line, then (len*n).p
		// must be between -tolerance and len+tolerance
		double xpr = dx*px + dy*py;
		if(xpr<0.) return false;
		if(xpr>distanceSq) return false;
	}
	
	else
	{	double dz = zmax-zmin, pz = v.z-zmin;
		
		// Let n = (dx,xy,dz)/len be unit vector along the line, then (len*n).p
		// must be between 0 and len^2
		double xpr = dx*px + dy*py + dz*pz;
		if(xpr<0.) return false;
		if(xpr>distanceSq) return false;
		
		// now get distance to the line d = |p - (p.n)n| = |p - xpr(dx,dy,dz)/distanceSq|
		xpr /= distanceSq;
		px -= xpr*dx;
		py -= xpr*dy;
		pz -= xpr*dz;
		double tpr = sqrt(px*px + py*py + pz*pz);
		if(tpr>tolerance) return false;
		if(tpr<-tolerance) return false;
	}
	
    return true;
}

// set a property
void LineController::SetProperty(const char *aName,char *value,CommonReadHandler *reader)
{
	if(strcmp(aName,"tolerance")==0)
	{	if(value[0]=='*')
		{	sscanf(&value[1],"%lf",&tolerance);
			tolerance*=ElementBase::GetMinimumCellSize();
		}
		else
		{	sscanf(value,"%lf",&tolerance);
			tolerance*=distScaling;
		}
	}
	else
		ShapeController::SetProperty(aName,value,reader);
}

#pragma mark LineController: Accessors

// set tolerance (no need of scaling)
void LineController::SetTolerance(double value) { tolerance=value; }

// type of object
const char *LineController::GetShapeName(void) { return "Line"; }


