/********************************************************************************
    OvalController.hpp
    nairn-mpm-fea

    Created by John Nairn on 4/23/12.
    Copyright (c) 2012 John A. Nairn, All rights reserved.

    Dependencies
        ShapeController.hpp
 ********************************************************************************/

#ifndef _OVALCONTROLLER_

#define _OVALCONTROLLER_

#include "Read_XML/ShapeController.hpp"

class OvalController : public ShapeController
{
    public:
	
        // contructors
        OvalController(int);
        virtual bool FinishSetup(void);
	
        // methods
        virtual bool ContainsPoint(Vector &);
    
        // accessors
        virtual const char *GetShapeName();
    
    private:
        double x0,y0;
};

#endif

