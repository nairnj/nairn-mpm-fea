/********************************************************************************
    BoxController.cpp
    NairnFEA and NairnMPM
    
    Created by John Nairn on 8/30/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_XML/BoxController.hpp"

/********************************************************************************
	LineController: Constructors and Destructor
********************************************************************************/

BoxController::BoxController(int block) : ShapeController(block)
{
}

/********************************************************************************
	LineController: methods
********************************************************************************/

// Deterime if point is sufficiently close to the line from
//   (xmin,ymin) to (xmax,ymax). Must be within rectangle a
//   distance tolerance from the line in all directions 
bool BoxController::PtOnShape(Vector v)
{	return v.x<=xmax && v.x>=xmin && v.y<=ymax && v.y>=ymin && v.z<=zmax && v.z>=zmin;
}
