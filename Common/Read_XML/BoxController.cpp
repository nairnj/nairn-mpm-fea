/********************************************************************************
    BoxController.cpp
    nairn-mpm-fea
    
    Created by John Nairn on 8/30/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_XML/BoxController.hpp"
#ifdef SUPPORT_MEMBRANES
    #include "MPM_Classes/MemPoint3D.hpp"
#endif

#pragma mark BoxController: constructors, destructor, initializers

BoxController::BoxController(int block) : ShapeController(block)
{
    axis = 0;			// 0 for box, 1,2,3 for axis of a cylinder
    coneRadius = 1.;    // between -1 and 1 for radius at top (>0) or botton (<0), cylinder only
#ifdef SUPPORT_MEMBRANES
	memaxis = 3;		// 1,2,3 for thickness direction of a membrane
#endif
	twoDShape = false;
}

// set a property
// throws SAXException()
void BoxController::SetProperty(const char *aName,char *value,CommonReadHandler *reader)
{
    if(strcmp(aName,"axis")==0)
    {	if(strcmp(value,"x")==0 || strcmp(value,"X")==0 || strcmp(value,"1")==0)
            axis=1;
        else if(strcmp(value,"y")==0 || strcmp(value,"Y")==0 || strcmp(value,"2")==0)
            axis=2;
        else if(strcmp(value,"z")==0 || strcmp(value,"Z")==0 || strcmp(value,"3")==0)
            axis=3;
        else if(strcmp(value,"0")==0)
            axis=0;
        else
            ThrowSAXException("Box axis must be x, y, z, 1, 2, or 3");
    }
    else if(strcmp(aName,"radius")==0)
    {   sscanf(value,"%lf",&coneRadius);
        if(coneRadius<-1. || coneRadius>1.)
            ThrowSAXException("Cone radius must be between -1 and 1");
    }
#ifdef SUPPORT_MEMBRANES
    // Membranes should add angle orientation in the future
    else if(strcmp(aName,"memaxis")==0)
    {	if(strcmp(value,"x")==0 || strcmp(value,"X")==0 || strcmp(value,"1")==0)
			memaxis=1;
		else if(strcmp(value,"y")==0 || strcmp(value,"Y")==0 || strcmp(value,"2")==0)
			memaxis=2;
        else if(strcmp(value,"z")==0 || strcmp(value,"Z")==0 || strcmp(value,"3")==0)
			memaxis=3;
		else
			ThrowSAXException("Box memaxis must be x, y, z, 1, 2, or 3");
    }
	else if(strcmp(aName,"length")==0)
	{	// membrane length
		sscanf(value,"%lf",&length);
		length*=distScaling;
	}
#endif // end SUPPORT_MEMBRANES
    else
        ShapeController::SetProperty(aName,value,reader);
}

// called after initialization is done
bool BoxController::FinishSetup(void)
{
    ShapeController::FinishSetup();
    
    // midpoints
    xmid = (xmin+xmax)/2.;
    ymid = (ymin+ymax)/2.;
    zmid = (zmin+zmax)/2.;
    
    // radii squared
    a2 = xmax-xmid;
    a2 *= a2;
    b2 = ymax-ymid;
    b2 *= b2;
    c2 = zmax-zmid;
    c2 *= c2;
	
    return TRUE;
}

#pragma mark BoxController: methods

// Deterine if point in box or cylinder
bool BoxController::ContainsPoint(Vector& v)
{
    if(axis==0)
    {   return v.x<=xmax && v.x>=xmin && v.y<=ymax && v.y>=ymin && v.z<=zmax && v.z>=zmin;
    }
    else if(axis==1)
    {   if(v.x>xmax || v.x<xmin) return FALSE;
        double dx = xmax-xmin;
        double R = coneRadius<0. ? (v.x-xmin+coneRadius*(v.x-xmax))/dx : (xmax-v.x+coneRadius*(v.x-xmin))/dx ;
        double dy = v.y-ymid;
        double dz = v.z-zmid;
        return (dy*dy/b2 + dz*dz/c2) <= R*R ;
    }
    else if(axis==2)
    {   if(v.y>ymax || v.y<ymin) return FALSE;
        double dy = ymax-ymin;
        double R = coneRadius<0. ? (v.y-ymin+coneRadius*(v.y-ymax))/dy : (ymax-v.y+coneRadius*(v.y-ymin))/dy ;
        double dx = v.x-xmid;
        double dz = v.z-zmid;
        return (dx*dx/a2 + dz*dz/c2) <= R*R ;
    }
    else
    {   if(v.z>zmax || v.z<zmin) return FALSE;
        double dz = zmax-zmin;
        double R = coneRadius<0. ? (v.z-zmin+coneRadius*(v.z-zmax))/dz : (zmax-v.z+coneRadius*(v.z-zmin))/dz ;
        double dx = v.x-xmid;
        double dy = v.y-ymid;
        return (dx*dx/a2 + dy*dy/b2) <= R*R ;
    }
}

#ifdef SUPPORT_MEMBRANES

// Membrane loop
#pragma mark BoxController: Membrane methods

// prepare to return particle locations for a membrane
bool BoxController::StartMembraneLoop(void)
{
	// false if membrane propertties not enetered
	if(length<=0. || a2<=0. || b2<=0. || c2<=0.) return false;
    coneRadius = 1.;                 // cones not allowed
	
	// x-y are virtial indicators for plane of the membrane
	double xsideLength,xRange,ysideLength,yRange;
	alpha = beta = gamma = 0.;
	if(memaxis==3)
	{	xsideLength = 2.*sqrt(a2);
		xRange = xmax-xmin;
		ysideLength = 2.*sqrt(b2);
		yRange = ymax-ymin;
		mthick = fabs(zmax-zmin);
	}
	else if(memaxis==2)
	{	xsideLength = 2.*sqrt(a2);
		xRange = xmax-xmin;
		ysideLength = 2.*sqrt(c2);
		yRange = zmax-zmin;
		mthick = fabs(ymax-ymin);
		alpha = beta = gamma = PI_CONSTANT/2.;
	}
	else
	{	xsideLength = 2.*sqrt(b2);
		xRange = ymax-ymin;
		ysideLength = 2.*sqrt(c2);
		yRange = zmax-zmin;
		mthick = fabs(xmax-xmin);
		beta = PI_CONSTANT/2.;
	}
	
	// get length per particle x side
	numXParticles = int(xsideLength/length+.5);
	if(numXParticles<1) numXParticles = 1;
	mxLength = xRange/(double)numXParticles;
	
    // y side
	numYParticles = int(ysideLength/length+.5);
	if(numYParticles<1) numYParticles = 1;
	myLength = yRange/(double)numYParticles;
    
    // total and start enumerator
    numParticles = numXParticles*numYParticles;
	particleNum = 0;
	return true;
}

// return next particle location
bool BoxController::NextMembraneLocation(Vector *loc)
{
	if(particleNum>=numParticles) return false;
    
    int row = int(particleNum/numXParticles);
    int col = particleNum % numXParticles;
    
	if(memaxis==3)
	{	loc->x = xmin + col*mxLength + mxLength/2.;
		loc->y = ymin + row*myLength + myLength/2.;
		loc->z = zmid;
	}
	else if(memaxis==2)
	{	loc->x = xmin + col*mxLength + mxLength/2.;
		loc->y = ymid;
		loc->z = zmin + row*myLength + myLength/2.;
	}
	else
	{	loc->x = xmid;
		loc->y = ymin + col*mxLength + mxLength/2.;
		loc->z = zmin + row*myLength + myLength/2.;
	}
    
    // in the future, rotate this point to final position
    
	particleNum++;
	return true;
}

// set properties on each material point
// z thickness, x length, and y length (membrance currently all in x-y plane
// with thickness direction in z direction - will change in the future).
void BoxController::SetMembraneProperties(MPMBase *mpt)
{
	((MemPoint3D *)mpt)->SetMembrane(mthick,mxLength,myLength,alpha,beta,gamma);
}

#endif // end SUPPORT_MEMBRANES

#pragma mark BoxController: accessors

// override for 3D objects
bool BoxController::Is2DShape(void) { return FALSE; }

// type of object
const char *BoxController::GetShapeName(void)
{   if(axis==0) return "Box";
    return coneRadius>-1. && coneRadius<1. ? "Cone" : "Cylinder" ;
}


