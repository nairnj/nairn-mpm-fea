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
	
		// methods
		virtual bool FinishSetup(void);
		virtual bool ContainsPoint(Vector &);
		virtual bool HasAllParameters(void);
		virtual void SetParameter(char *,char *);
		virtual void FinishParameter(void);
	
	private:
		double xparm,yparm;
		vector< double > xpt;
		vector< double > ypt;
};

#endif


