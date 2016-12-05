/********************************************************************************
    EdgeBCController.cpp
    NairnFEA
    
    Created by John Nairn on 6/25/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_FEA/EdgeBCController.hpp"
#include "Elements/ElementBase.hpp"
#include "Boundary_Conditions/EdgeBC.hpp"

EdgeBCController *edgeBCCtrl=NULL;

/********************************************************************************
	EdgeBCController: methods
********************************************************************************/

// add object to linked list, but only if unique
bool EdgeBCController::AddObject(LinkedObject *obj)
{
	// If already set, delete and return FALSE
	EdgeBC *newObject=(EdgeBC *)obj;
	EdgeBC *nextObject=(EdgeBC *)firstObject;
	while(nextObject!=NULL)
	{	if(nextObject->SameDofSetting(newObject))
		{	delete newObject;
			return FALSE;
		}
		nextObject=(EdgeBC *)nextObject->GetNextObject();
	}
	
	// add new one
	return ParseController::AddObject(obj);
}

// set stress values for current edge BC
int EdgeBCController::SetStress(char *xData)
{
	double stress[3];
	
	if(lastObject==NULL) return FALSE;
	
	int elemID=theElements[((EdgeBC *)lastObject)->ElementIndex()]->FaceNodes();
	switch(elemID)
	{	case 2:
			sscanf(xData,"%lf,%lf",&stress[0],&stress[1]);
			break;
		case 3:
			 sscanf(xData,"%lf,%lf,%lf",&stress[0],&stress[1],&stress[2]);
			break;
	   default:
			return FALSE;
	}
	((EdgeBC *)lastObject)->SetStress(stress,elemID);
	return TRUE;
}
