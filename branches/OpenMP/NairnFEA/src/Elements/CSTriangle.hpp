/********************************************************************************
    CSTriangle.hpp
    NairnFEA
    
    Created by John Nairn on 10/24/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
	
	Dependencies
		Linear2D.hpp (ElementBase.hpp)
********************************************************************************/

#ifndef _CSTRIANGLE_

#define _CSTRIANGLE_

#include "Elements/Linear2D.hpp"

class CSTriangle : public Linear2D
{
    public:
        // constructors
        CSTriangle(int,int *,int,double,double);
        
        // prototypes
        virtual short ElementName(void);
		void Stiffness(int);
		void ForceStress(double *,int,int);
	
		// const methods
		virtual int NumberNodes(void) const;
		virtual int NumberSides(void) const;
		virtual void ShapeFunction(Vector *,int,double *,double *,double *,
									Vector *,double *,double *,double *) const;
};

#endif
