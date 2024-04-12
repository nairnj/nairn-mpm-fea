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
#include "Materials/MaterialBase.hpp"

class MatPtFluxBC : public MatPtLoadBC
{
    public:
		int face;
        int phaseStyle;
    
        // constructors and destructors
		MatPtFluxBC(int,int,int,int,int);
		
		// methods
		virtual BoundaryCondition *PrintBC(ostream &);
        MatPtLoadBC *AddMPFluxBC(double);
	
		// accessors
		virtual void SetBCValue(double);
		virtual void SetBCFirstTime(double);
		virtual double GetBCFirstTimeOut(void);
};

extern MatPtFluxBC *firstFluxPt;
extern MatPtFluxBC *firstDiffFluxBC[NUM_DUFFUSION_OPTIONS];

#endif
