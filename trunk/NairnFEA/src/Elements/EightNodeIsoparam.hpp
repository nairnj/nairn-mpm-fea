/********************************************************************************
    EightNodeIsoparam.hpp
    NairnFEA
    
    Created by John Nairn on Fri Oct 22 2002.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
	
	Dependencies
		Quad2D.hpp (ElementBase.hpp)
********************************************************************************/

#ifndef _EIGHTNODEISOPARAM_

#define _EIGHTNODEISOPARAM_

#include "Elements/Quad2D.hpp"

class EightNodeIsoparam : public Quad2D
{
    public:
        // constructors
        EightNodeIsoparam(int,int *,int,double,double);
        
        // prototypes
        virtual short ElementName(void);
        virtual int NumberNodes(void);
        virtual void ShapeFunction(Vector *,int,double *,double *,double *,
                                Vector *,double *,double *,double *);
};

#endif

