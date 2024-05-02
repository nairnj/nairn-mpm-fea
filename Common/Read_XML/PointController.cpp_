/********************************************************************************
    PointController.cpp
    NairnFEA
    
    Created by John Nairn on 8/8/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Read_XML/PointController.hpp"

/********************************************************************************
	LineController: Constructors and Destructor
********************************************************************************/

PointController::PointController(int block,int node) : ShapeController(block)
{
	sourceBlock=block;
	nearestNode=node;
}

// bothing to check
bool PointController::FinishSetup(void) { return TRUE; }

/********************************************************************************
	LineController: methods
********************************************************************************/

// return next node for this line
int PointController::nextNode(void)
{
	if(nodeNum>1) return 0;
	nodeNum++;
	return nearestNode;
}

#ifdef MPM_CODE
// not used for particles
int PointController::nextParticle(void) { return -1; }
#endif

// type of object
const char *PointController::GetShapeName(void) { return "Point"; }
