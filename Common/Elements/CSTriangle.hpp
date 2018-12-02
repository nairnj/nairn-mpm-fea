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
#ifdef MPM_CODE
		CSTriangle(int,int *);
#else
        CSTriangle(int,int *,int,double,double);
#endif
        
        // prototypes
        virtual short ElementName(void);
		virtual void ShapeFunction(Vector *,int,double *,double *,double *,
								   Vector *,double *,double *,double *) const;
#ifdef MPM_CODE
		virtual void ShapeFunction(Vector *,int,double *,double *,double *,double *) const;

#else
		void Stiffness(int);
		void ForceStress(double *,int,int);
#endif
	
		// const methods
		virtual int NumberNodes(void) const;
		virtual int NumberSides(void) const;
};

#endif
