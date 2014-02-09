/********************************************************************************
    LinkedObject.hpp
    NairnFEA
    
    Created by John Nairn on 6/22/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _LINKEDOBJECT_

#define _LINKEDOBJECT_

class LinkedObject
{
    public:
		LinkedObject *nextObject;
	
        //  Constructors and Destructor
        LinkedObject();
		
		// methods
		void SetNextObject(LinkedObject *);
		LinkedObject *GetNextObject(void) const;
};

#endif
