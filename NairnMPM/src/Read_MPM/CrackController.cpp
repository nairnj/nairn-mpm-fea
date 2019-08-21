/********************************************************************************
    CrackController.cpp
    NairnFEA
    
    Created by John Nairn on 6/27/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Read_MPM/CrackController.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Cracks/CrackSegment.hpp"

CrackController *crackCtrl=NULL;

#pragma mark CrackController: methods

// start new crack
void CrackController::AddCrack(CrackHeader *newCrack)
{
	int newNumber = (lastObject==NULL) ? 1 :
					((CrackHeader *)lastObject)->GetNumber()+1;
	newCrack->SetNumber(newNumber);
	AddObject(newCrack);
}

// Add new segment
int CrackController::AddSegment(CrackSegment *newSegment,bool ignoreDuplicate)
{	
	return ((CrackHeader *)lastObject)->addSegment(newSegment,ignoreDuplicate);
}

// When done with crack, finish any other needed initializations
bool CrackController::FinishCrack(void)
{
	return ((CrackHeader *)lastObject)->CreateHierarchy();
}

// assemble into array used in the code
int CrackController::SetCracksArray(void)
{
	// make 0-based array of cracks
	if(numObjects==0) return false;
	crackList = new (std::nothrow) CrackHeader *[numObjects];
	if(crackList==NULL) return false;
	
	// fill the array
	CrackHeader *obj = (CrackHeader *)firstObject;
	numberOfCracks = 0;
	while(obj!=NULL)
	{	crackList[numberOfCracks]=obj;
		numberOfCracks++;
		obj=(CrackHeader *)obj->GetNextObject();
	}
 	return true;
}


