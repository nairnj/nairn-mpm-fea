/********************************************************************************
    LinkedObject.hpp
    NairnFEA
    
    Created by John Nairn on 6/22/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"

/********************************************************************************
	LinkedObject: Constructors and Destructor
********************************************************************************/

LinkedObject::LinkedObject()
{
	nextObject=NULL;
}

/********************************************************************************
	LinkedObject: methods
********************************************************************************/

// next object
LinkedObject *LinkedObject::GetNextObject(void) const { return nextObject; }
void LinkedObject::SetNextObject(LinkedObject *obj) { nextObject=obj; }
