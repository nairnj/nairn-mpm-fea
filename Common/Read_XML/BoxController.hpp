/********************************************************************************
    BoxController.hpp
    NairnFEA and NairnMPM
    
    Created by John Nairn on 8/30/07.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
	
	Dependencies
		ShapeController.hpp
********************************************************************************/

#ifndef _BOXCONTROLLER_

#define _BOXCONTROLLER_

#include "Read_XML/ShapeController.hpp"

class BoxController : public ShapeController
{
    public:
        int axis;
	
		// contructors
		BoxController(int);
	
		// methods
        virtual void SetProperty(const char *,char *,CommonReadHandler *);
        virtual void FinishSetup(void);
		virtual bool PtOnShape(Vector);
    
    protected:
        double xmid,ymid,zmid;
        double a2,b2,c2;
		
};

#endif
