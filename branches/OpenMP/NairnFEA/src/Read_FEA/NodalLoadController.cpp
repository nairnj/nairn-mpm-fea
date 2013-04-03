/********************************************************************************
    NodalLoadController.cpp
    NairnFEA
    
    Created by John Nairn on 10/18/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_FEA/NodalLoadController.hpp"
#include "Boundary_Conditions/NodalLoad.hpp"

NodalLoadController *loadBCCtrl=NULL;

/********************************************************************************
	NodalDispBCController: methods
********************************************************************************/

// add object to linked list, but only if unique
bool NodalLoadController::AddObject(LinkedObject *obj)
{
	// If already set, delete and return FALSE
	NodalLoad *newObject=(NodalLoad *)obj;
	NodalLoad *nextObject=(NodalLoad *)firstObject;
	while(nextObject!=NULL)
	{	if(nextObject->SameDofSetting(newObject))
		{	delete newObject;
			return FALSE;
		}
		nextObject=(NodalLoad *)nextObject->GetNextObject();
	}
	
	// add new one
	return ParseController::AddObject(obj);
}

