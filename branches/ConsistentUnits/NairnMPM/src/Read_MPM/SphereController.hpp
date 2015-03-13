/********************************************************************************
    SphereController.hpp
    NairnFEA
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		ShapeController.hpp
********************************************************************************/

#ifndef _SPHERECONTROLLER_

#define _SPHERECONTROLLER_

#include "Read_XML/ShapeController.hpp"

class SphereController : public ShapeController
{
    public:
	
        // initialize
        SphereController(int);
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



