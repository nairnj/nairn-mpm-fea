/******************************************************************
	EightNodeIsoparamBrick.hpp
	nairn-mpm-fea

	Created by John Nairn on 7/20/06.
	Copyright 2006 RSAC Software. All rights reserved.
	
	Dependencies
		ElementBase3D.hpp (ElementBase.hpp)
 ******************************************************************/

#ifndef _EIGHTNODEISOPARAMBRICK_

#define _EIGHTNODEISOPARAMBRICK_

#include "Elements/ElementBase3D.hpp"

class EightNodeIsoparamBrick : public ElementBase3D
{
    public:
        // constructors
		EightNodeIsoparamBrick(int,int *);
        
        // prototypes
        virtual short ElementName(void);
		short PtInElement(Vector &);
		virtual int Orthogonal(double *,double *,double *);
	
		// const methods
		virtual int NumberNodes(void) const;
		virtual double GetArea(void) const;
		virtual double GetVolume(void) const;
		virtual double GetThickness(void) const;
	
		virtual void GetGimpNodes(int *,int *,int *,Vector *,Vector &) const;
		virtual void GimpShapeFunction(Vector *,int,int *,int,double *,double *,double *,double *,Vector &) const;
		virtual void GetXiPos(Vector *,Vector *) const;
		virtual void ShapeFunction(Vector *,int,double *,double *,double *,
										Vector *,double *,double *,double *) const;
		virtual void ShapeFunction(Vector *,int,double *,double *,double *,double *) const;
};

#endif
