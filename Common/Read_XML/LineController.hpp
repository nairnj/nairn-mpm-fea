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
#ifdef MPM_CODE
	class MPMBase;
#endif

class LineController : public ShapeController
{
    public:
	
		// contructors
		LineController(int,bool);
		LineController(int,double,double,double,double,double);
	
        // initialize
        virtual void SetProperty(const char *,char *,CommonReadHandler *);
        virtual bool FinishSetup(void);
    
		// methods
		virtual bool ContainsPoint(Vector &);
		
		// line only methods
		void SetTolerance(double);

#ifdef SUPPORT_MEMBRANES
		virtual bool StartMembraneLoop(void);
		virtual bool NextMembraneLocation(Vector *);
		virtual void SetMembraneProperties(MPMBase *);
#endif
    
        // accessors
        virtual const char *GetShapeName();

	protected:
		double distanceSq,dtolerance,tolerance;
#ifdef SUPPORT_MEMBRANES
		double thickness,length,sineAngle;
		double xLength,yLength;
#endif
	
};

#endif
