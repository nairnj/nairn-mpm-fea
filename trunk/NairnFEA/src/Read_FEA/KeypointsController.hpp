/********************************************************************************
    KeypointsController.hpp
    NairnFEA
    
    Created by John Nairn on 6/22/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		ParceController.hpp
********************************************************************************/

#ifndef _KEYPOINTS_

#define _KEYPOINTS_

#include "Read_XML/ParseController.hpp"

class Keypoint;

class KeypointsController : public ParseController
{
    public:
	
		~KeypointsController(void);
		
		// methods
		int ValidName(char *);
		Keypoint *FindKeypoint(char *);

};

extern KeypointsController *keyPts;

#endif
