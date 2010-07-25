/********************************************************************************
    BodyBoxController.hpp
    NairnFEA
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		BodyObjectController.hpp
********************************************************************************/

#ifndef _BODYBOXCONTROLLER_

#define _BODYBOXCONTROLLER_

#include "Read_MPM/BodyObjectController.hpp"

class BodyBoxController : public BodyObjectController
{
    public:
	
		// methods
		virtual bool ContainsPoint(Vector &);
		virtual bool Is2DBodyObject(void);
};

#endif


