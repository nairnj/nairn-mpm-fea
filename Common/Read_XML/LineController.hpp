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
	
		// methods
		virtual bool PtOnShape(Vector);
		virtual void SetProperty(char *,char *,CommonReadHandler *);
		virtual void FinishSetup(void);
		
		// line only methods
		void SetTolerance(double);

	protected:
		double distanceSq,dtolerance,tolerance;
	
};

#endif
