/********************************************************************************
    NodalDispBCController.hpp
    NairnFEA
    
    Created by John Nairn on 10/18/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		ParceController.hpp
********************************************************************************/

#ifndef _NODALDISPBCCONTROLLER_

#define _NODALDISPBCCONTROLLER_

#include "Read_XML/ParseController.hpp"

class NodalDispBCController : public ParseController
{
    public:
	
		// methods
		bool AddObject(LinkedObject *);

};

extern NodalDispBCController *dispBCCtrl;

#endif
