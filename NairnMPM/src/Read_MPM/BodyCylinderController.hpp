/********************************************************************************
    BodyCylinderController.hpp
    NairnFEA
    
    Created by John Nairn on 9/1/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		BodyObjectController.hpp
		BodySphereController.hpp
********************************************************************************/

#ifndef _BODYCYLINDERCONTROLLER_

#define _BODYCYLINDERCONTROLLER_

#include "Read_MPM/BodySphereController.hpp"

class BodyCylinderController : public BodySphereController
{
    public:
	
		// constructor
		BodyCylinderController();
		
		// methods
		virtual bool FinishSetup(void);
		virtual bool ContainsPoint(Vector &);
		virtual void SetBodyProperty(const char *,char *,CommonReadHandler *);
	
	protected:
		int axis;			// 1, 2, or 3, for x, y or z
	
};

#endif

