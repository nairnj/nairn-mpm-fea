/********************************************************************************
    QuadInterface.hpp
    NairnFEA
    
    Created by John Nairn on 01/07/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
	
	Dependencies
		Interface2D.hpp (ElementBase.hpp)
********************************************************************************/

#ifndef _QUADINTERFACE_

#define _QUADINTERFACE_

#include "Elements/Interface2D.hpp"

class QuadInterface : public Interface2D
{
    public:
        // constructors
		QuadInterface(long,long *,int,double,double); 
        
        // prototypes
        virtual short ElementName(void);
        virtual int NumberNodes(void);
        virtual void ShapeFunction(Vector *,int,double *,double *,double *,
                                    Vector *,double *,double *,double *);
		void Stiffness(int);
		virtual int FaceNodes(void);
};

#endif
