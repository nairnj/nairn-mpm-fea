/********************************************************************************
    TorusController.cpp
    nairn-mpm-fea

    Created by John Nairn on 8/30/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_MPM/TorusController.hpp"

#pragma mark TorusController: constructors, destructor, initializers

TorusController::TorusController(int block) : ShapeController(block)
{
    axis = 3;           // 0 for box, 1,2,3 for axis of a torus (normal to plane)
    ringRadius = -1.;   // if not set use axis radius
}

// set a property
// throws SAXException()
void TorusController::SetProperty(const char *aName,char *value,CommonReadHandler *reader)
{
    if(strcmp(aName,"axis")==0)
    {	if(strcmp(value,"x")==0 || strcmp(value,"X")==0 || strcmp(value,"1")==0)
            axis=1;
        else if(strcmp(value,"y")==0 || strcmp(value,"Y")==0 || strcmp(value,"2")==0)
            axis=2;
        else if(strcmp(value,"z")==0 || strcmp(value,"Z")==0 || strcmp(value,"3")==0)
            axis=3;
        else
            ThrowSAXException("Torus axis must be x, y, z, 1, 2, or 3");
    }
    else if(strcmp(aName,"radius")==0)
    {   sscanf(value,"%lf",&ringRadius);
        if(ringRadius<=0.)
            ThrowSAXException("Torus ring radius must be positive");
    }
    else
        ShapeController::SetProperty(aName,value,reader);
}

// called after initialization is done
bool TorusController::FinishSetup(void)
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
    
    if(ringRadius>0.)
        d2 = ringRadius*ringRadius;
    else if(axis==1)
        d2 = a2;
    else if(axis==2)
        d2 = b2;
    else
        d2 = c2;
    
    return TRUE;
}

#pragma mark TorusController: methods

// Deterine if point in box or cylinder
bool TorusController::ContainsPoint(Vector& v)
{
	// R is distance along line from plane midpoint to plane elipse
	// rpos is distance from plane midpoint to the point in the plane
	double R;
	
    if(axis==1)
    {   if(v.x>xmax || v.x<xmin) return FALSE;
		if(DbleEqual(v.y-ymid,0.))
			R = sqrt(c2);
		else
		{	double m = (v.z-zmid)/(v.y-ymid);
			R = sqrt((1.+m*m)/(1./b2 + m*m/c2));
		}
        double rpos = sqrt((v.z-zmid)*(v.z-zmid)+(v.y-ymid)*(v.y-ymid));
        return (R-rpos)*(R-rpos)/d2 + (v.x-xmid)*(v.x-xmid)/a2 <= 1.;
    }
    else if(axis==2)
    {   if(v.y>ymax || v.y<ymin) return FALSE;
		if(DbleEqual(v.x-xmid,0.))
			R = sqrt(c2);
		else
		{	double m = (v.z-zmid)/(v.x-xmid);
			R = sqrt((1.+m*m)/(1./a2 + m*m/c2));
		}
        double rpos = sqrt((v.x-xmid)*(v.x-xmid)+(v.z-zmid)*(v.z-zmid));
        return (R-rpos)*(R-rpos)/d2 + (v.y-ymid)*(v.y-ymid)/b2 <= 1.;
    }
    else
    {   if(v.z>zmax || v.z<zmin) return FALSE;
		if(DbleEqual(v.x-xmid,0.))
			R = sqrt(b2);
		else
		{	double m = (v.y-ymid)/(v.x-xmid);
			R = sqrt((1.+m*m)/(1./a2 + m*m/b2));
		}
        double rpos = sqrt((v.x-xmid)*(v.x-xmid)+(v.y-ymid)*(v.y-ymid));
        return (R-rpos)*(R-rpos)/d2 + (v.z-zmid)*(v.z-zmid)/c2 <= 1.;
    }
}

#pragma mark TorusController: accessors

// override for 3D objects
bool TorusController::Is2DShape(void) { return FALSE; }

// type of object
const char *TorusController::GetShapeName(void) { return "Torus"; }

