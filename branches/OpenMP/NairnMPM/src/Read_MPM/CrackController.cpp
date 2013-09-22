/********************************************************************************
    CrackController.cpp
    NairnFEA
    
    Created by John Nairn on 6/27/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_MPM/CrackController.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Cracks/CrackSegment.hpp"

CrackController *crackCtrl=NULL;

/********************************************************************************
	CrackController: methods
********************************************************************************/

// start new crack
void CrackController::AddCrack(CrackHeader *newCrack)
{
	int newNumber = (lastObject==NULL) ? 1 :
					((CrackHeader *)lastObject)->GetNumber()+1;
	newCrack->SetNumber(newNumber);
	AddObject(newCrack);
}

// Add new segment
int CrackController::AddSegment(CrackSegment *newSegment)
{
	return ((CrackHeader *)lastObject)->add(newSegment);
}

// When done with crack, finish any other needed initializations
void CrackController::FinishCrack(void)
{
#ifdef HIERARCHICAL_CRACKS
	((CrackHeader *)lastObject)->CreateHierarchy();
#endif
}


