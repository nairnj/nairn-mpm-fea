/********************************************************************************
    LineController.cpp
    NairnFEA
    
    Created by John Nairn on 6/62/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Elements/ElementBase.hpp"
#include "Read_XML/LineController.hpp"
#ifdef SUPPORT_MEMBRANES
	#include "MPM_Classes/MemPoint2D.hpp"
	#include "MPM_Classes/MemPointAS.hpp"
#endif

#pragma mark LineController: Constructors and Destructor

// Create, may be 2D or 3D line
LineController::LineController(int block,bool is2D) : ShapeController(block)
{
	tolerance = 0.;
#ifdef SUPPORT_MEMBRANES
	thickness = -1;
	length = -1.;
#endif
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
#ifdef SUPPORT_MEMBRANES
	else if(strcmp(aName,"thickness")==0)
	{	// membrane thickness
		sscanf(value,"%lf",&thickness);
		thickness*=distScaling;
	}
	else if(strcmp(aName,"length")==0)
	{	// membrane length
		sscanf(value,"%lf",&length);
		length*=distScaling;
	}
#endif
	else
		ShapeController::SetProperty(aName,value,reader);
}

#ifdef SUPPORT_MEMBRANES

#pragma mark LineController: Membrane methods

// prepare to return particle locations for a membrane
bool LineController::StartMembraneLoop(void)
{
	// false if membrane propertties not enetered
	if(length<=0. || thickness<=0. || distanceSq==0.) return false;
	
	// get length per particle
	double lineLength = sqrt(distanceSq);
	numParticles = int(lineLength/length+.5);
	if(numParticles<1) numParticles = 1;
	length = lineLength/(double)numParticles;
	
	// start enuerator
	particleNum = 0;
	xLength = (xmax-xmin)/(double)numParticles;
	yLength = (ymax-ymin)/(double)numParticles;
	
	// sine of CW rotation angle of x axis, but between -90 (to +y) and 90 (to -y)
	// Rotation matrix = {{cos(angle), sineAngle},{-sineAngle,cos(angle)}
	if(xLength>0)
		sineAngle = (ymin-ymax)/lineLength;
	else
		sineAngle = (ymax-ymin)/lineLength;
	
	return true;
}

// return next particle location
bool LineController::NextMembraneLocation(Vector *loc)
{
	if(particleNum>=numParticles) return false;
	loc->x = xmin + particleNum*xLength + xLength/2.;
	loc->y = ymin + particleNum*yLength + yLength/2.;
	particleNum++;
	return true;
}

// set properties on each material point
// length, thickness, and sine(cw angle)
void LineController::SetMembraneProperties(MPMBase *mpt)
{
	// this will call MemPointAS if needed
	((MemPoint2D *)mpt)->SetMembrane(thickness,length,sineAngle);
}

#endif // end SUPPORT_MEMBRANES

#pragma mark LineController: Accessors

// set tolerance (no need of scaling)
void LineController::SetTolerance(double value) { tolerance=value; }

// type of object
const char *LineController::GetShapeName(void) { return "Line"; }


