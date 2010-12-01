/******************************************************************
	EightNodeIsoparamBrick.hpp
	NairnMPM

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
        virtual int NumberNodes(void);
		virtual double GetArea(void);
		virtual double GetVolume(void);
		virtual double GetThickness(void);
		short PtInElement(Vector &);
        virtual void ShapeFunction(Vector *,int,double *,double *,double *,
                                Vector *,double *,double *,double *);
		virtual void ShapeFunction(Vector *,int,double *,double *,double *,double *);
		virtual int Orthogonal(double *,double *,double *);
		virtual void GetGimpNodes(int *,int *,int *,Vector *);
		virtual void GimpShapeFunction(Vector *,int,int *,int,double *,double *,double *,double *);
		virtual void GetXiPos(Vector *,Vector *);
};

#endif
