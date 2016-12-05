/********************************************************************************
    RectController.hpp
    NairnFEA
    
    Created by John Nairn on 8/8/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		ShapeController.hpp
********************************************************************************/

#ifndef _RECTCONTROLLER_

#define _RECTCONTROLLER_

#include "Read_XML/ShapeController.hpp"

class RectController : public ShapeController
{
    public:
	
		// contructors
		RectController(int);
		virtual void SetParameter(const char *,const char *);
		virtual bool FinishParameter(void);
	
		// methods
		virtual bool ContainsPoint(Vector &);
		virtual bool CheckAngle(double,double);
	
        // accessors
        virtual const char *GetShapeName();
	
	protected:
		double startAngle;
		double endAngle;


};

#endif

