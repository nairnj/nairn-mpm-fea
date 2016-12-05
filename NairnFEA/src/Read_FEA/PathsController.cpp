/********************************************************************************
    PathsController.cpp
    NairnFEA
    
    Created by John Nairn on 6/23/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_FEA/PathsController.hpp"
#include "Read_FEA/Path.hpp"

PathsController *paths=NULL;

/********************************************************************************
	PathsController: contructor and destructor
********************************************************************************/

// delete keypoints
PathsController::~PathsController(void)
{
	// delete the linked paths
	Path *prevPath,*aPath=(Path *)firstObject;
	while(aPath!=NULL)
	{	prevPath=aPath;
		aPath=(Path *)prevPath->GetNextObject();
		delete prevPath;
	}
}

/********************************************************************************
	PathsController: methods
********************************************************************************/

// add key point to the last path
int PathsController::AddKeypoint(char *keyName)
{
	if(lastObject==NULL) return FALSE;
	return ((Path *)lastObject)->AddKeypoint(keyName);
}

// check if last path now valid
int PathsController::ValidPath(void)
{
	if(lastObject==NULL) return FALSE;
	return ((Path *)lastObject)->ValidPath();
}

// Check name before use
int PathsController::ValidName(char *pathName)
{
	if(pathName==NULL) return FALSE;
	int len=(int)strlen(pathName);
	if(len<=0 || len>32) return FALSE;
	if(FindPath(pathName)!=NULL) return FALSE;
	return TRUE;
}

// get pointer to a path
Path *PathsController::FindPath(char *pathName)
{
	if(pathName[0]==0) return NULL;
	Path *path=(Path *)firstObject;
	while(path!=NULL)
	{	if(strcmp(pathName,path->pathID)==0)
			return path;
		path=(Path *)path->GetNextObject();
	}
	return NULL;
}

