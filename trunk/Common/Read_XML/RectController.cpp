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

/********************************************************************************
	RectController: methods
********************************************************************************/

// Deterime if point is sufficiently close to the line from
//   (xmin,ymin) to (xmax,ymax). Must be within rectangle a
//   distance tolerance from the line in all directions 
bool RectController::ContainsPoint(Vector& v)
{	return v.x<=xmax && v.x>=xmin && v.y<=ymax && v.y>=ymin;
}

// type of object
const char *RectController::GetShapeName(void) { return "Rect"; }
