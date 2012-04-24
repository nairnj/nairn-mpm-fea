/********************************************************************************
    LineController.hpp
    NairnFEA
    
    Created by John Nairn on 6/26/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		ShapeController.hpp
********************************************************************************/

#ifndef _LINECONTROLLER_

#define _LINECONTROLLER_

#include "Read_XML/ShapeController.hpp"

class CommonReadHandler;

class LineController : public ShapeController
{
    public:
	
		// contructors
		LineController(int);
		LineController(int,double,double,double,double,double);
	
        // initialize
        virtual void SetProperty(const char *,char *,CommonReadHandler *);
        virtual bool FinishSetup(void);
    
		// methods
		virtual bool ContainsPoint(Vector &);
		
		// line only methods
		void SetTolerance(double);
    
        // accessors
        virtual const char *GetShapeName();

	protected:
		double distanceSq,dtolerance,tolerance;
	
};

#endif
