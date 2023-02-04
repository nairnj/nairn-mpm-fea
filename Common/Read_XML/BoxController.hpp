/********************************************************************************
    BoxController.hpp
    nairn-mpm-fea
    
    Created by John Nairn on 8/30/07.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
	
	Dependencies
		ShapeController.hpp
********************************************************************************/

#ifndef _BOXCONTROLLER_

#define _BOXCONTROLLER_

#include "Read_XML/ShapeController.hpp"

#ifdef MPM_CODE
    class MPMBase;
#endif

class BoxController : public ShapeController
{
    public:
        int axis;
        double coneRadius;
	
		// contructors
		BoxController(int);
	
        // initialize
        virtual void SetProperty(const char *,char *,CommonReadHandler *);
        virtual bool FinishSetup(void);
    
		// methods
		virtual bool ContainsPoint(Vector &);

#ifdef SUPPORT_MEMBRANES
        virtual bool StartMembraneLoop(void);
        virtual bool NextMembraneLocation(Vector *);
        virtual void SetMembraneProperties(MPMBase *);
#endif
        // accessors
        virtual bool Is2DShape(void);
        virtual const char *GetShapeName();
    
    protected:
        double xmid,ymid,zmid;
        double a2,b2,c2;
#ifdef SUPPORT_MEMBRANES
		int memaxis;
        double length,mxLength,myLength,mthick;
        int numXParticles,numYParticles;
		double alpha,beta,gamma;
#endif
		
};

#endif
