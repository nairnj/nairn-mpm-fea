/********************************************************************************
    ParseController.cpp
    NairnFEA
    
    Created by John Nairn on 6/22/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_XML/ParseController.hpp"

#pragma mark ParseController: Constructors and Destructor

ParseController::ParseController()
{
	firstObject=NULL;
	lastObject=NULL;
	numObjects=0;
}

#pragma mark ParseController: methods

// add object to linked list
bool ParseController::AddObject(LinkedObject *obj)
{
	numObjects++;
	if(firstObject==NULL)
		firstObject=(LinkedObject *)obj;
	else
		lastObject->SetNextObject((LinkedObject *)obj);
	lastObject=(LinkedObject *)obj;
	return TRUE;
}

// last object is also th e current object when reading XMP file
LinkedObject *ParseController::currentObject(void) { return lastObject; }

// delete all linked objects
void ParseController::ClearObjects(void)
{
	while(firstObject!=NULL)
	{	LinkedObject *prevObject = firstObject;
		firstObject = prevObject->GetNextObject();
		delete prevObject;
	}
	firstObject=NULL;
	lastObject=NULL;
	numObjects=0;
}


