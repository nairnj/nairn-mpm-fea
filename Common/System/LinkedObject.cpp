/********************************************************************************
    LinkedObject.hpp
    NairnFEA
    
    Created by John Nairn on 6/22/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"

#pragma mark LinkedObect::Constructors and Destructor

LinkedObject::LinkedObject()
{
	nextObject=NULL;
}

#pragma mark LinkedObect::Methods

// next object
LinkedObject *LinkedObject::GetNextObject(void) const { return nextObject; }
void LinkedObject::SetNextObject(LinkedObject *obj) { nextObject=obj; }

// hack to get previous object when it is not stored
// return NULL if is firstObject or -1L if is not in the list or if firstObject is NULL
LinkedObject *LinkedObject::GetPreviousObject(LinkedObject *firstObject) const
{
	// if this is firstObject (which cannot be NULL) return NULL
	if(this==firstObject) return NULL;
	
	// scan until nextObject is this
	while(firstObject!=NULL)
	{	LinkedObject *cmpObject = firstObject->GetNextObject();
		if(this==cmpObject) return firstObject;
		firstObject = cmpObject;
	}
	
	// getting here means this is not in the list
	return (LinkedObject *)-1L;
}

#pragma mark LinkedObject::Class Methods

// Delete deleteObj in the list starting with *firstObj
// Return object to one after the one that was deleted (or NULL if was the last one)
// Not very efficient; if deleting many objects, better to write custom code
// If not in the list (or firstObject is NULL), nothing is changed or deleted, but still returns the next object
// Not currently used. Better test first
LinkedObject *LinkedObject::DeleteObject(LinkedObject *deleteObj,LinkedObject **firstObj)
{
	// get next and previous one (note that deleteObj cannot be NULL)
	LinkedObject *nextObj = deleteObj->GetNextObject();
	LinkedObject *prevObj = deleteObj->GetPreviousObject(*firstObj);
	
	// if not found, do nothing, but return next one anyway
	if(prevObj==(LinkedObject *)-1L)
		return nextObj;
	
	// if was first, change first, otherwise set next of previous objects
	if(prevObj==NULL)
		*firstObj = nextObj;
	else
		prevObj->SetNextObject(nextObj);
	
	// delete and return next one
	delete deleteObj;
	return nextObj;
}

// insert addObj object after afterObj and return the added one for convenience
// if firstObj!=NULL, both addObj and afterObj must be valid pointers
// Odd behavior - addObj is change on call to next memory location. Not sure why currently not user
LinkedObject *LinkedObject::InsertObject(LinkedObject *addObj,LinkedObject *afterObj,LinkedObject **firstObj)
{
	if(*firstObj==NULL)
	{	// if empty, insert at start of linked list
		*firstObj = addObj;
	}
	else
	{	// set after to point to added on and added on to point to one after previous afterOBj
		LinkedObject *oldNextObj = afterObj->GetNextObject();		// NULL if was last object
		afterObj->SetNextObject(addObj);
		addObj->SetNextObject(oldNextObj);							// NULL if was last
	}
	return addObj;
}


