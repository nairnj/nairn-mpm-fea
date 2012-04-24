/********************************************************************************
    BodyRectController.hpp
    NairnFEA
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		BodyObjectController.hpp
********************************************************************************/

#ifndef _BODYRECTCONTROLLER_

#define _BODYRECTCONTROLLER_

#include "Read_MPM/BodyObjectController.hpp"

class BodyRectController : public BodyObjectController
{
    public:
	
		// methods
		virtual bool ContainsPoint(Vector &);
    
        // accessors
		virtual const char *GetShapeName(void);
};

#endif

