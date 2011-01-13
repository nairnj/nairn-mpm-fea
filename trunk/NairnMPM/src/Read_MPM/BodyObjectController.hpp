/********************************************************************************
    BodyObjectController.hpp
    NairnFEA
    
    Created by John Nairn on 8/10/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _BODYOBJECTCONTROLLER_

#define _BODYOBJECTCONTROLLER_

class CommonReadHandler;

class BodyObjectController
{
    public:
	
		// contructors
		BodyObjectController();
		virtual ~BodyObjectController();
	
		// methods
		virtual bool ContainsPoint(Vector &);
		virtual bool FinishSetup(void);
		virtual bool Is2DBodyObject(void);
		virtual void SetBodyProperty(const char *,char *,CommonReadHandler *);
		virtual void SetParameter(char *,char *);
		virtual bool FinishParameter(void);
		virtual void SetScaling(double);
		virtual bool HasAllParameters(void);
		virtual bool SetBodyPropertyFromData(char *,CommonReadHandler *);
	
	protected:
		double xmin,xmax,ymin,ymax,zmin,zmax;
		double distScaling;
};

extern BodyObjectController *theBody;

#endif
