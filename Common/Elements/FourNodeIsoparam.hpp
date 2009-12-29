/********************************************************************************
    FourNodeIsoparam.hpp
    NairnMPM
    
    Created by John Nairn on Wed Jan 24 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		Linear2D.hpp (ElementBase.hpp)
********************************************************************************/

#ifndef _FOURNODEISOPARAM_

#define _FOURNODEISOPARAM_

#include "Elements/Linear2D.hpp"

class FourNodeIsoparam : public Linear2D
{
    public:
        // constructors
#ifdef MPM_CODE
        FourNodeIsoparam(long,long *);
#else
        FourNodeIsoparam(long,long *,int,double,double);
#endif
        
        // prototypes
        virtual short ElementName(void);
        virtual int NumberNodes(void);
        virtual void ShapeFunction(Vector *,int,double *,double *,double *,
                                Vector *,double *,double *,double *);
#ifdef MPM_CODE
		virtual void ShapeFunction(Vector *,int,double *,double *,double *,double *);
        virtual void FindExtent(void);
		virtual int Orthogonal(double *,double *,double *);
		virtual void GetGimpNodes(int *,int *,int *,Vector *);
		virtual void GimpShapeFunction(Vector *,int,int *,int,double *,double *,double *,double *);
#endif

	protected:
#ifdef MPM_CODE
        virtual void GetXiPos(Vector *,Vector *);
#endif

};

#endif

