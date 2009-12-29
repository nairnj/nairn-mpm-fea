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
		virtual bool PtOnShape(Vector);
		virtual void SetProperty(char *,char *,CommonReadHandler *);
		virtual void FinishSetup(void);
	
	private:
		double startAngle,endAngle;
		
	protected:
		double centerX, centerY, a, b;

};

#endif

