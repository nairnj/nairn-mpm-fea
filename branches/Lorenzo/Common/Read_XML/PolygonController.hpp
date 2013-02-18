/********************************************************************************
    PolygonController.hpp
    NairnFEA
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		ShapeController.hpp
********************************************************************************/

#ifndef _POLYGONCONTROLLER_

#define _POLYGONCONTROLLER_

#include "Read_XML/ShapeController.hpp"

class PolygonController : public ShapeController
{
    public:
	
		// constructors
        PolygonController(int);
		virtual ~PolygonController();
	
        // initialize
        virtual bool FinishSetup(void);
        virtual void SetParameter(const char *,const char *);
        virtual bool FinishParameter(void);
        virtual bool HasAllParameters(void);
    
		// methods
		virtual bool ContainsPoint(Vector &);
    
        // accessors
		virtual const char *GetShapeName(void);
	
	private:
		double xparm,yparm;
		vector< double > xpt;
		vector< double > ypt;
};

#endif


