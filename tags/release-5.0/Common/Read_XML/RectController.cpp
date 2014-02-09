/********************************************************************************
    RectController.cpp
    NairnFEA
    
    Created by John Nairn on 8/8/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_XML/RectController.hpp"

#pragma mark RectController: constructors, destructor, initializers

RectController::RectController(int block) : ShapeController(block)
{
}

#pragma mark RectController: methods

// return true if point is in this rectangle
bool RectController::ContainsPoint(Vector& v)
{	return v.x<=xmax && v.x>=xmin && v.y<=ymax && v.y>=ymin;
}

#pragma mark RectController: accessors

// type of object
const char *RectController::GetShapeName(void) { return "Rectangle"; }
