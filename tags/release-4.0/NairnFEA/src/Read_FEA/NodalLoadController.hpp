/********************************************************************************
    NodalLoadController.hpp
    NairnFEA
    
    Created by John Nairn on 10/18/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		ParceController.hpp
********************************************************************************/

#ifndef _NODALLOADCONTROLLER_

#define _NODALLOADCONTROLLER_

#include "Read_XML/ParseController.hpp"

class NodalLoadController : public ParseController
{
    public:
	
		// methods
		bool AddObject(LinkedObject *);

};

extern NodalLoadController *loadBCCtrl;

#endif
