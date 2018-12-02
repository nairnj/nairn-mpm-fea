/********************************************************************************
	Lagrange2D.hpp
	NairnFEA

	Created by John Nairn on 7 Feb 2011.
	Copyright (c) 2011 John A. Nairn, All rights reserved.

	Dependencies
		Quad2D.hpp (ElementBase.hpp)
********************************************************************************/

#ifndef _LAGRANGE2D_

#define _LAGRANGE2D_

#include "Elements/Quad2D.hpp"

class Lagrange2D : public Quad2D
{
    public:
        // constructors
#ifdef MPM_CODE
		Lagrange2D(int,int *);
#else
		Lagrange2D(int,int *,int,double,double);
#endif
	
		// prototypes
		virtual short ElementName(void);
		virtual int NumberNodes(void) const;
		virtual double GetArea(void) const;
		virtual void ShapeFunction(Vector *,int,double *,double *,double *,
							   Vector *,double *,double *,double *) const;
#ifdef FEA_CODE
		virtual void ExtrapolateGaussStressToNodes(double [][5]);
#endif
#ifdef MPM_CODE
		virtual void ShapeFunction(Vector *,int,double *,double *,double *,double *) const;
		virtual void SplineShapeFunction(int *,Vector *,int,double *,double *,double *,double *) const;
		virtual void FindExtent(void);
		virtual int Orthogonal(double *,double *,double *);
	
		// const methods
		virtual void GetXiPos(const Vector *,Vector *) const;
#endif
	
};

#endif
		
