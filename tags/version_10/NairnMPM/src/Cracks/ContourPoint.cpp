/********************************************************************************
    ContourPoint.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Feb 13 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Cracks/ContourPoint.hpp"
#include "Nodes/NodalPoint.hpp"

/*******************************************************************
	ContourPoint: Constructors and Destructors
*******************************************************************/

// Constructors
ContourPoint::ContourPoint()
{
}

ContourPoint::ContourPoint(NodalPoint *aNode)
{
    node=aNode;
    nextPoint=NULL;
    orient=ANGLED;
	phantomNode = false;
}

ContourPoint::~ContourPoint()
{
	// delete the node if it was a phantom node
	if(phantomNode)
	{	delete node;
	}
}

/*******************************************************************
	ContourPoint: Methods
*******************************************************************/

/* set next point and find orientation
    return ANGLED (-1) if line not along x or y axis
*/
int ContourPoint::SetNextPoint(ContourPoint *apt)
{
    double dl;
    
    nextPoint=apt;
    double dx=nextPoint->node->x-node->x;
    double dy=nextPoint->node->y-node->y;
    if(DbleEqual(dx,0.))
    {	orient=VERTICAL;
        norm.x=dy>0 ? 1. : -1.;
        norm.y=0;
        ds=dy/norm.x;
    }
    else if(DbleEqual(dy,0.))
    {	orient=HORIZONTAL;
        norm.x=0.;
        norm.y=dx>0 ? -1. : 1.;
        ds=-dx/norm.y;
    }
    else
    {	orient=ANGLED;
        dl=sqrt(dx*dx+dy*dy);
        norm.x=dy/dl;
        norm.y=-dx/dl;
        ds=dl;			// redo if use off-axis contours
    }
    return orient;
}

// find relative location on this segment (assume on the segment to nextPoint)
double ContourPoint::Fraction(Vector &pt)
{
    double fraction=0.;
    
    switch(orient)
    {	case HORIZONTAL:
            fraction=(pt.x-node->x)/(nextPoint->node->x-node->x);
            break;
        
        case VERTICAL:
            fraction=(pt.y-node->y)/(nextPoint->node->y-node->y);
            break;
    }
    return fraction;
}

// mark as a phantom node that was inserted at point where crack crosses the contour
void ContourPoint::SetPhantomNode(bool phantom) { phantomNode = phantom; }

