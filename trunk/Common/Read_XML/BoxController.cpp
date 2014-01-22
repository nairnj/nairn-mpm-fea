/********************************************************************************
    BoxController.cpp
    nairn-mpm-fea
    
    Created by John Nairn on 8/30/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_XML/BoxController.hpp"

#pragma mark BoxController: constructors, destructor, initializers

BoxController::BoxController(int block) : ShapeController(block)
{
    axis = 0;    // 0 for box, 1,2,3 for axis of a cylinder
    coneRadius = 1.;    // between -1 and 1 for radius at top (>0) or botton (<0), cylinder only
}

// set a property
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
    else
        ShapeController::SetProperty(aName,value,reader);
}

// called after initialization is done
bool BoxController::FinishSetup(void)
{
    ShapeController::FinishSetup();
    
    // done for box
    if(axis==0) return TRUE;
    
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

#pragma mark BoxController: accessors

// override for 3D objects
bool BoxController::Is2DShape(void) { return FALSE; }

// type of object
const char *BoxController::GetShapeName(void)
{   if(axis==0) return "Box";
    return coneRadius>-1. && coneRadius<1. ? "Cone" : "Cylinder" ;
}


