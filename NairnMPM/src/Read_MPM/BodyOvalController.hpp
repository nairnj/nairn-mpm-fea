/********************************************************************************
    BodyOvalController.hpp
    NairnFEA
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		BodyObjectController.hpp
********************************************************************************/

#ifndef _BODYOVALCONTROLLER_

#define _BODYOVALCONTROLLER_

#include "Read_MPM/BodyObjectController.hpp"

class BodyOvalController : public BodyObjectController
{
    public:
	
		// methods
		virtual bool FinishSetup(void);
		virtual bool ContainsPoint(Vector &);
		virtual char *GetObjectType(void);
	
	private:
		double x0,y0;
};

#endif


