/********************************************************************************
    PointController.hpp
    NairnFEA
    
    Created by John Nairn on 8/8/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		ShapeController.hpp
********************************************************************************/

#ifndef _POINTCONTROLLER_

#define _POINTCONTROLLER_

#include "Read_XML/ShapeController.hpp"

class PointController : public ShapeController
{
    public:
	
		// contructors
		PointController(int,int);
        virtual bool FinishSetup(void);
	
		// methods
		virtual int nextNode(void);
#ifdef MPM_CODE
		virtual int nextParticle(void);
#endif
    
        // accessors
        virtual const char *GetShapeName();

	private:
		int nearestNode;
	
};

#endif
