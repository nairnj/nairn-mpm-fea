/********************************************************************************
    KeypointsController.cpp
    NairnFEA
    
    Created by John Nairn on 6/22/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Read_FEA/KeypointsController.hpp"
#include "Read_FEA/Keypoint.hpp"

KeypointsController *keyPts=NULL;

/********************************************************************************
	KeypointsController: contructor and destructor
********************************************************************************/

// delete keypoints
KeypointsController::~KeypointsController(void)
{
	// delete the linked keypoints
	Keypoint *prevKey,*aKey=(Keypoint *)firstObject;
	while(aKey!=NULL)
	{	prevKey=aKey;
		aKey=(Keypoint *)prevKey->GetNextObject();
		delete prevKey;
	}
}

/********************************************************************************
	KeypointsController: methods
********************************************************************************/

// Check name before use
int KeypointsController::ValidName(char *keyName)
{
	if(keyName==NULL) return FALSE;
	int len=strlen(keyName);
	if(len<=0 || len>32) return FALSE;
	if(FindKeypoint(keyName)!=NULL) return FALSE;
	return TRUE;
}

// get pointer to a keypoint
Keypoint *KeypointsController::FindKeypoint(char *keyName)
{
	if(keyName[0]==0) return NULL;
	Keypoint *key=(Keypoint *)firstObject;
	while(key!=NULL)
	{	if(strcmp(keyName,key->keyID)==0)
			return key;
		key=(Keypoint *)key->GetNextObject();
	}
	return NULL;
}

