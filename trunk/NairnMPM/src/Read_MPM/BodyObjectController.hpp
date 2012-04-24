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
	
        // Initialize
        virtual void SetProperty(const char *,char *,CommonReadHandler *);
        virtual void SetProperty(const char *,double);
        virtual void SetProperty(char *,CommonReadHandler *);
        virtual void SetParameter(const char *,const char *);
        virtual bool FinishParameter(void);
        virtual bool HasAllParameters(void);
        virtual void SetScaling(double);
        virtual bool FinishSetup(void);
    
		// methods
		virtual bool ContainsPoint(Vector &);
    
        // accessors
		virtual bool Is2DShape(void);
		virtual const char *GetShapeName();
	
	protected:
		double xmin,xmax,ymin,ymax,zmin,zmax;
		double distScaling;
};

extern BodyObjectController *theBody;

#endif
