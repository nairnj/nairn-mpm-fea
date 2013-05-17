/********************************************************************************
    TorusController.hpp
    NairnMPM

    Created by John Nairn on 5/15/13.
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Dependencies
        ShapeController.hpp
********************************************************************************/

#ifndef _TORUSCONTROLLER_

#define _TORUSCONTROLLER_

#include "Read_XML/ShapeController.hpp"

class TorusController : public ShapeController
{
    public:
        int axis;
        double ringRadius;
        
        // contructors
        TorusController(int);
        
        // initialize
        virtual void SetProperty(const char *,char *,CommonReadHandler *);
        virtual bool FinishSetup(void);
        
        // methods
        virtual bool ContainsPoint(Vector &);
        
        // accessors
        virtual bool Is2DShape(void);
        virtual const char *GetShapeName();
        
    protected:
        double xmid,ymid,zmid;
        double a2,b2,c2,d2;
        
};

#endif
