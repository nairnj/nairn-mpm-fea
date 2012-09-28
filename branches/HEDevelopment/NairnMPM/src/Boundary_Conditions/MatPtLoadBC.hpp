/********************************************************************************
    MatPtLoadBC.hpp
    NairnMPM
    
    Created by John Nairn on Tues Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		BoundaryCondition.hpp
********************************************************************************/

#ifndef _MATPTLOADBC_

#define _MATPTLOADBC_

#include "Boundary_Conditions/BoundaryCondition.hpp"

class MatPtLoadBC: public BoundaryCondition
{
    public:
        int ptNum;
        int direction;
        
        // constructors and destructors
        MatPtLoadBC(int,int,int);
        
        // virtual methods
        virtual BoundaryCondition *PrintBC(ostream &);
		virtual void GetPosition(double *,double *,double *,double *);
		
		// specific methods
        MatPtLoadBC *ZeroMPLoad(void);
        MatPtLoadBC *AddMPLoad(double);
        MatPtLoadBC *ReverseLinearLoad(double,double *);
        MatPtLoadBC *MakeConstantLoad(double);
        
		// class methods
		static void SetParticleFext(double);
};

extern MatPtLoadBC *firstLoadedPt;

#endif
