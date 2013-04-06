/********************************************************************************
    LinearInterface.hpp
    NairnFEA
    
    Created by John Nairn on 01/07/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
	
	Dependencies
		Interface2D.hpp (ElementBase.hpp)
********************************************************************************/

#ifndef _LINEARINTERFACE_

#define _LINEARINTERFACE_

#include "Elements/Interface2D.hpp"

class LinearInterface : public Interface2D
{
    public:
        // constructors
		LinearInterface(int,int *,int,double,double); 
        
        // prototypes
        virtual short ElementName(void);
        virtual int NumberNodes(void) const;
        virtual void ShapeFunction(Vector *,int,double *,double *,double *,
                                    Vector *,double *,double *,double *) const;
		virtual int FaceNodes(void);
};

#endif
