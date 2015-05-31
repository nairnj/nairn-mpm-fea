/********************************************************************************
    MatPtFluxBC.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Mar 17 2004.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		MatPtLoadBC.hpp (BoundaryCondition.hpp)
********************************************************************************/

#ifndef _MATPTFLUXBC_

#define _MATPTFLUXBC_

#include "Boundary_Conditions/MatPtLoadBC.hpp"

class MatPtFluxBC : public MatPtLoadBC
{
    public:
		int face;
    
        // constructors and destructors
		MatPtFluxBC(int,int,int,int);
		
		// methods
		BoundaryCondition *PrintBC(ostream &);
        MatPtFluxBC *AddMPFlux(double);
		virtual void SetBCValue(double);
		virtual void SetBCFirstTime(double);
		virtual double GetBCFirstTimeOut(void);
};

extern MatPtFluxBC *firstFluxPt;

#endif
