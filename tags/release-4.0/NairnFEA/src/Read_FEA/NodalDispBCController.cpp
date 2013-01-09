/********************************************************************************
    NodalDispBCController.cpp
    NairnFEA
    
    Created by John Nairn on 10/18/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_FEA/NodalDispBCController.hpp"
#include "Boundary_Conditions/NodalDispBC.hpp"

NodalDispBCController *dispBCCtrl=NULL;

/********************************************************************************
	NodalDispBCController: methods
********************************************************************************/

// add object to linked list, but only if unique
bool NodalDispBCController::AddObject(LinkedObject *obj)
{
	// If already set, delete and return FALSE
	NodalDispBC *newObject=(NodalDispBC *)obj;
	NodalDispBC *nextObject=(NodalDispBC *)firstObject;
	while(nextObject!=NULL)
	{	if(nextObject->SameDofSetting(newObject))
		{	delete newObject;
			return FALSE;
		}
		nextObject=(NodalDispBC *)nextObject->GetNextObject();
	}
	
	// add new one
	return ParseController::AddObject(obj);
}

