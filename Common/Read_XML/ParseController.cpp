/********************************************************************************
    ParseController.cpp
    NairnFEA
    
    Created by John Nairn on 6/22/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Read_XML/ParseController.hpp"

#pragma mark ParseController: Constructors and Destructor

#include "Common/System/LinkedObject.hpp"
#include "System/MPMPrefix.hpp"

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

// delete all linked objects starting with one object
void ParseController::ClearObjects(LinkedObject *obj)
{
    // if first one, clear them all
    if(obj==firstObject)
    {   ClearObjects();
        return;
    }
    
    bool deleting = false;
    numObjects = 0;
    LinkedObject *nextObj = firstObject;
    while(nextObj!=NULL)
    {   LinkedObject *holdObj = nextObj->GetNextObject();
        if(deleting)
        {   delete nextObj;
        }
        else if(holdObj==obj)
        {   // now reached the starting object to start deleting
            nextObj->SetNextObject(NULL);
            lastObject = nextObj;
            deleting = true;
        }
        else
        {   // count those remaining
            numObjects++;
        }
        nextObj = holdObj;
    }
}
