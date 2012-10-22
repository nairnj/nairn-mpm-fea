/********************************************************************************
    ArcController.hpp
    NairnFEA
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		LineController.hpp,ShapeController.hpp
********************************************************************************/

#ifndef _ARCCONTROLLER_

#define _ARCCONTROLLER_

#include "Read_XML/LineController.hpp"

class CommonReadHandler;

class ArcController : public LineController
{
    public:
	
		// contructors
		ArcController(int);
	
		// methods
		virtual bool ContainsPoint(Vector &);
		virtual void SetProperty(const char *,char *,CommonReadHandler *);
		virtual bool FinishSetup(void);
    
        // accessors
        virtual const char *GetShapeName();
	
	private:
		double startAngle,endAngle;
		
	protected:
		double centerX, centerY, a, b;

};

#endif

