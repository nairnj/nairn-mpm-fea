/********************************************************************************
    SixNodeTriangle.hpp
    NairnFEA
    
    Created by John Nairn on Mon Oct 25 2002.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
	
	Dependencies
		Quad2D.hpp (ElementBase.hpp)
********************************************************************************/

#ifndef _SIXNODETRIANGLE_

#define _SIXNODETRIANGLE_

#include "Elements/Quad2D.hpp"

class SixNodeTriangle : public Quad2D
{
    public:
        // constructors
        SixNodeTriangle(int,int *,int,double,double);
        
        // prototypes
        virtual short ElementName(void);
        virtual void ShapeFunction(Vector *,int,double *,double *,double *,
                                Vector *,double *,double *,double *) const;
		virtual void ExtrapolateGaussStressToNodes(double [][5]);
	
		// const methods
		virtual int NumberNodes(void) const;
		virtual int NumberSides(void) const;
};

#endif

