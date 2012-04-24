/********************************************************************************
    BodyPolygonController.hpp
    NairnFEA
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		BodyObjectController.hpp
********************************************************************************/

#ifndef _BODYPOLYGONCONTROLLER_

#define _BODYPOLYGONCONTROLLER_

#include "Read_MPM/BodyObjectController.hpp"

class BodyPolygonController : public BodyObjectController
{
    public:
	
		// constructors
		~BodyPolygonController();
	
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


