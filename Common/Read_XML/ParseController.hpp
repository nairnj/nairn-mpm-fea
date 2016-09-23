/********************************************************************************
    ParseController.hpp
    NairnFEA
    
    Created by John Nairn on 6/22/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _PARSECONTROLLER_

#define _PARSECONTROLLER_

class ParseController
{
    public:
		LinkedObject *firstObject,*lastObject;
		int numObjects;
	
        //  Constructors and Destructor
        ParseController();
		
		// methods
		bool AddObject(LinkedObject *);
		LinkedObject *currentObject(void);
};

#endif
