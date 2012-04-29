/********************************************************************************
    Linear2D.hpp
    NairnMPM
    
    Created by John Nairn on Wed Jan 24 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		ElementBase.hpp
********************************************************************************/

#ifndef _LINEAR2D_

#define _LINEAR2D_

#include "Elements/ElementBase.hpp"

class Linear2D : public ElementBase
{
    public:
        double thickness;
        
        // constructors
#ifdef MPM_CODE
        Linear2D(int,int *);
#else
        Linear2D(int,int *,int,double,double);
#endif
        
        // prototypes
		virtual int FaceNodes(void);
        virtual double GetArea(void);
		virtual double GetVolume(void);
        virtual short PtInElement(Vector &);
        virtual double GetThickness(void);
    virtual void SetThickness(double);
#ifdef FEA_CODE
		virtual void CalcEdgeLoads(double *,int,int,double *,int);
#endif
};

#endif
