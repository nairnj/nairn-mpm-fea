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
    
        // initialize
        virtual void SetProperty(const char *,char *,CommonReadHandler *);
        virtual void SetProperty(const char *,double);
        virtual void SetProperty(char *,CommonReadHandler *);
        virtual void SetScaling(double);
        virtual void SetParameter(const char *,const char *);
        virtual bool FinishParameter(void);
        virtual bool FinishSetup(void);
        virtual bool HasAllParameters(void);
    
		// methods
		virtual bool ContainsPoint(Vector &);
		virtual void resetNodeEnumerator(void);
		virtual const char *startNodeEnumerator(int,int);
		virtual int nextNode(void);
        void resetElementEnumerator(void);
        int nextElement(void);
    
        // MPM only methods
#ifdef MPM_CODE
		void setNetBC(bool);
		double particleCount(void);
		virtual void resetParticleEnumerator(void);
		virtual int nextParticle(void);
#endif
    
        // accessors
        virtual const char *GetShapeName(void);
        virtual bool Is2DShape(void);
        virtual char *GetContextInfo(void);
		
		virtual double Returnx0(void); //modiftf ******** #rigidbodyrotation
		virtual double Returny0(void); //modiftf ******** #rigidbodyrotation
		virtual double Returnz0(void); //modiftf ******** #rigidbodyrotation
		
		
        // base class only (non virtual)
        int GetSourceBlock(void);
        bool RequiredBlock(int);
		
	protected:
		double xmin,xmax,ymin,ymax,zmin,zmax;
		double distScaling;
		int sourceBlock,nodeNum,elemNum;
#ifdef MPM_CODE
		int particleNum,numParticles;
#endif
		
};

extern ShapeController *theShape;

#endif

