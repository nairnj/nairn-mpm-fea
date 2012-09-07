/********************************************************************************
    PathBCController.cpp
    NairnFEA
    
    Created by John Nairn on  8/8/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_FEA/PathBCController.hpp"
#include "Read_FEA/Path.hpp"

#pragma mark PathBCController: initilizers

PathBCController::PathBCController(int block,Path *bcPath) : ShapeController(block)
{
	thePath=bcPath;
}

// bothing to check
bool PathBCController::FinishSetup(void) { return TRUE; }

#pragma mark PathBCController: methods

// reset nodeNum and verify OK if applied to a path in FEA
const char *PathBCController::startNodeEnumerator(int command,int axis)
{
	nodeNum=1;
	return thePath->VerifyBCOption(command,axis);
}

// return next node for this line
int PathBCController::nextNode(void)
{
	nodeNum++;
	return thePath->nodeAtIndex(nodeNum-1);
}

#pragma mark PathBCController: accessors

// return the path
char *PathBCController::GetContextInfo(void) { return (char *)thePath; }

// type of object
const char *PathBCController::GetShapeName(void) { return "Path"; }
