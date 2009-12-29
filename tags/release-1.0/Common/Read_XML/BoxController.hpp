/********************************************************************************
    BoxController.hpp
    NairnFEA and NairnMPM
    
    Created by John Nairn on 8/30/07.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
	
	Dependencies
		ShapeController.hpp
********************************************************************************/

#ifndef _BOXCONTROLLER_

#define _BOXCONTROLLER_

#include "Read_XML/ShapeController.hpp"

class BoxController : public ShapeController
{
    public:
	
		// contructors
		BoxController(int);
	
		// methods
		virtual bool PtOnShape(Vector);
		
};

#endif
