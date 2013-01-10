/********************************************************************************
    FEAArchiveData.hpp
    NairnFEA
    
    Created by John Nairn on Jan 15, 2006.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
	
	Dependencies
		CommonArchiveData.hpp
********************************************************************************/

#ifndef _FEAARCHIVEDATA_

#define _FEAARCHIVEDATA_

#include "System/CommonArchiveData.hpp"

class FEAArchiveData : public CommonArchiveData
{
    public:
	
        //  Constructors and Destructor
        FEAArchiveData();
		
		// abstract methods
	
		// methods

		// accessors
};

extern FEAArchiveData *archiver;

#endif
