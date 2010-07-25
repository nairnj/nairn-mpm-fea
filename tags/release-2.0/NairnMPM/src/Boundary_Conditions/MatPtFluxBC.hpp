/********************************************************************************
    MatPtFluxBC.hpp
    NairnMPM
    
    Created by John Nairn on Wed Mar 17 2004.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		MatPtLoadBC.hpp (BoundaryCondition.hpp)
********************************************************************************/

#ifndef _MATPTFLUXDBC_

#define _MATPTFLUXBC_

#include "Boundary_Conditions/MatPtLoadBC.hpp"

class MatPtFluxBC : public MatPtLoadBC
{
    public:
        
        // constructors and destructors
		MatPtFluxBC(long,int,int);
		
		// methods
		BoundaryCondition *PrintBC(ostream &);
		MatPtFluxBC *ZeroMPFlux(void);
        MatPtFluxBC *AddMPFlux(double);
};

extern MatPtFluxBC *firstFluxPt;

#endif
