/********************************************************************************
    BoxController.cpp
    NairnFEA and NairnMPM
    
    Created by John Nairn on 8/30/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_XML/BoxController.hpp"

/********************************************************************************
	BoxController: Constructors and Destructor
********************************************************************************/

BoxController::BoxController(int block) : ShapeController(block)
{
    axis = 0;    // 0 for box, 1,2,3 for axis of a cylinder
}

/********************************************************************************
	BoxController: methods
********************************************************************************/

// Deterime if point is sufficiently close to the line from
//   (xmin,ymin) to (xmax,ymax). Must be within rectangle a
//   distance tolerance from the line in all directions 
bool BoxController::ContainsPoint(Vector& v)
{   if(axis==0)
    {   return v.x<=xmax && v.x>=xmin && v.y<=ymax && v.y>=ymin && v.z<=zmax && v.z>=zmin;
    }
    else if(axis==1)
    {   if(v.x>xmax || v.x<xmin) return FALSE;
        double dy = v.y-ymid;
        double dz = v.z-zmid;
        return (dy*dy/b2 + dz*dz/c2) <= 1. ;
    }
    else if(axis==2)
    {   if(v.y>ymax || v.y<ymin) return FALSE;
        double dx = v.x-xmid;
        double dz = v.z-zmid;
        return (dx*dx/a2 + dz*dz/c2) <= 1. ;
    }
    else
    {   if(v.z>zmax || v.z<zmin) return FALSE;
        double dx = v.x-xmid;
        double dy = v.y-ymid;
        return (dx*dx/a2 + dy*dy/b2) <= 1. ;
    }
}
// set a property
void BoxController::SetProperty(const char *aName,char *value,CommonReadHandler *reader)
{	if(strcmp(aName,"axis")==0)
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
    
    // radii
    a2 = xmax-xmid;
    a2 *= a2;
    b2 = ymax-ymid;
    b2 *= b2;
    c2 = zmax-zmid;
    c2 *= c2;
    
    return TRUE;
}

// override for 3D objects
bool BoxController::Is2DShape(void) { return FALSE; }

// type of object
const char *BoxController::GetShapeName(void) { return axis==0 ? "Box" : "Cylinder" ; }


