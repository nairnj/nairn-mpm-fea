/********************************************************************************
    BodySphereController.hpp
    NairnFEA
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		BodyObjectController.hpp
********************************************************************************/

#ifndef _BODYSPHERECONTROLLER_

#define _BODYSPHERECONTROLLER_

#include "Read_MPM/BodyObjectController.hpp"

class BodySphereController : public BodyObjectController
{
    public:
	
        // initialize
        virtual bool FinishSetup(void);
    
		// methods
		virtual bool ContainsPoint(Vector &);
    
        // accessors
		virtual bool Is2DShape(void);
		virtual const char *GetShapeName(void);
	
	protected:
		double x0,y0,z0;
};

#endif



