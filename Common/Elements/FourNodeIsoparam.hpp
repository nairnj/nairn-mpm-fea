/********************************************************************************
    FourNodeIsoparam.hpp
    nairn-mpm-fea
    
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
        FourNodeIsoparam(int,int *);
#else
        FourNodeIsoparam(int,int *,int,double,double);
#endif
        
        // prototypes
        virtual short ElementName(void);
        virtual int NumberNodes(void) const;
        virtual void ShapeFunction(Vector *,int,double *,double *,double *,
                                Vector *,double *,double *,double *) const;
#ifdef FEA_CODE
		virtual void ExtrapolateGaussStressToNodes(double [][5]);
#endif
#ifdef MPM_CODE
        virtual void FindExtent(void);
		virtual int Orthogonal(double *,double *,double *);
		virtual void GetPosition(Vector *xipos,Vector *);
		virtual short PtInElement(Vector &) const;
	
		// const methods
        virtual void ShapeFunction(Vector *,int,double *,double *,double *,double *) const;
		virtual void SplineShapeFunction(int *,Vector *,int,double *,double *,double *,double *) const;
		virtual void GimpShapeFunction(Vector *,int *,int,double *,double *,double *,double *,Vector &) const;
        virtual void TartanGimpShapeFunction(Vector *,int *,int,double *,double *,double *,double *,Vector &) const;
		virtual void GimpShapeFunctionAS(Vector *,int *,int,double *,double *,double *,double *,Vector &) const;
		virtual void BGimpShapeFunction(Vector *,int *,int,double *,double *,double *,double *,Vector &) const;
		virtual void BGimpShapeFunctionAS(Vector *,int *,int,double *,double *,double *,double *,Vector &) const;
		virtual void GetXiPos(const Vector *,Vector *) const;
		virtual void GridTractionFunction(Vector *,Vector *,bool,double *,int *,double,double) const;
		virtual void FiniteGimpShapeFunction(int *,double *,double *,double *,double *,MPMBase *)  const;
#endif

};

#endif

