/********************************************************************************
    RectController.cpp
    NairnFEA
    
    Created by John Nairn on 8/8/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_XML/RectController.hpp"

/********************************************************************************
	RectController: Constructors and Destructor
********************************************************************************/

RectController::RectController(int block) : ShapeController(block)
{
}

// called after initialization is done
void RectController::FinishSetup(void)
{
    double temp;
    if(xmin>xmax)
	{	temp=xmax;
        xmax=xmin;
        xmin=temp;
    }
    if(DbleEqual(xmin,xmax))
        ThrowSAXException("xmax can not equal xmin in input parameters.");
		
    if(ymin>ymax)
	{	temp=ymax;
        ymax=ymin;
        ymin=temp;
    }
    if(DbleEqual(ymin,ymax))
        ThrowSAXException("ymax can not equal ymin in input parameters.");
}

/********************************************************************************
	RectController: methods
********************************************************************************/

// Deterime if point is sufficiently close to the line from
//   (xmin,ymin) to (xmax,ymax). Must be within rectangle a
//   distance tolerance from the line in all directions 
bool RectController::PtOnShape(Vector v)
{	return v.x<=xmax && v.x>=xmin && v.y<=ymax && v.y>=ymin;
}

