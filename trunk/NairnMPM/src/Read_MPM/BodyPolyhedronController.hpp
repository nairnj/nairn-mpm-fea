/********************************************************************************
	BodyPolyhedronController.hpp
	NairnMPM

	Created by John Nairn on 1/7/11.
	Copyright (c) 2011 John A. Nairn, All rights reserved.

	Dependencies
		BodyObjectController.hpp
********************************************************************************/

#ifndef _BODYPOLYHEDRONCONTROLLER_

#define _BODYPOLYHEDRONCONTROLLER_

#include "Read_MPM/BodyObjectController.hpp"

class PolyTriangle;

enum { NO_FACES=0,TRICLINIC_POINTS,TRICLINIC_VECTORS,BOX_CORNERS } ;

class BodyPolyhedronController : public BodyObjectController
{
	public:
		// constructors
		~BodyPolyhedronController();
	
		// methods
		virtual bool FinishSetup(void);
		virtual bool ContainsPoint(Vector &);
		virtual bool HasAllParameters(void);
		virtual void SetParameter(const char *,const char *);
		virtual bool FinishParameter(void);
		virtual bool Is2DBodyObject(void);
		virtual bool SetBodyPropertyFromData(char *,CommonReadHandler *);
		virtual const char *GetObjectType(void);
	
	private:
		int style;
		vector< PolyTriangle * > faces;
		Vector pmin,pmax;
		char order[9];
	
};

#endif

