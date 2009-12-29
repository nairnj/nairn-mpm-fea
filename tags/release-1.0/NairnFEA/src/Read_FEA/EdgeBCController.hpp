/********************************************************************************
    EdgeController.hpp
    NairnFEA
    
    Created by John Nairn on 6/25/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		ParceController.hpp
********************************************************************************/

#ifndef _EDGEBCCONTROLLER_

#define _EDGEBCCONTROLLER_

#include "Read_XML/ParseController.hpp"

class EdgeBCController : public ParseController
{
    public:
	
		// methods
		bool AddObject(LinkedObject *obj);
		int SetStress(char *xData);

};

extern EdgeBCController *edgeBCCtrl;

#endif
