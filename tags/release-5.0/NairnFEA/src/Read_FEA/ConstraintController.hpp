/********************************************************************************
    ConstrainController.hpp
    NairnFEA
    
    Created by John Nairn on 2/6/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		ParceController.hpp
********************************************************************************/

#ifndef _CONSTRAINTCONTROLLER_

#define _CONSTRAINTCONTROLLER_

#include "Read_XML/ParseController.hpp"

class ConstraintController : public ParseController
{
    public:
	
		// methods
		bool AddObject(LinkedObject *);

};

extern ConstraintController *constraintCtrl;

#endif
