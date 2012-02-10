/********************************************************************************
    ShapeController.hpp
    NairnFEA
    
    Created by John Nairn on 8/8/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _SHAPECONTROLLER_

#define _SHAPECONTROLLER_

class CommonReadHandler;

class ShapeController
{
	public:
	
		// contructors
		ShapeController(int);
		ShapeController(int,double,double,double,double,double);
		virtual ~ShapeController();
    
        // base close only (non virtual)
        int GetSourceBlock(void);
        bool RequiredBlock(int);
		
		// methods
		virtual bool PtOnShape(Vector);
		virtual void resetNodeEnumerator(void);
		virtual const char *startNodeEnumerator(int,int);
		virtual int nextNode(void);
		virtual char *GetContextInfo(void);
		virtual void SetScaling(double);
		virtual void SetProperty(const char *,char *,CommonReadHandler *);
		virtual void SetProperty(const char *,double);
		virtual void FinishSetup(void);
#ifdef MPM_CODE
		void setNetBC(bool);
		double particleCount(void);
		virtual void resetParticleEnumerator(void);
		virtual int nextParticle(void);
#endif
		
	protected:
		double xmin,xmax,ymin,ymax,zmin,zmax;
		double distScaling;
		int sourceBlock,nodeNum;
#ifdef MPM_CODE
		int particleNum,numParticles;
#endif
		
};

extern ShapeController *theShape;

#endif

