/********************************************************************************
    MatPtheatFluxBC.hpp
    nairn-mpm-fea

    Created by John Nairn on May 28, 2013.
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Dependencies
        MatPtLoadBC.hpp (BoundaryCondition.hpp)
********************************************************************************/

#ifndef _MATPTHEATFLUXDBC_

#define _MATPTHEATFLUXDBC_

#include "Boundary_Conditions/MatPtLoadBC.hpp"

class MatPtHeatFluxBC : public MatPtLoadBC
{
    public:
        int face;
    
        // constructors and destructors
        MatPtHeatFluxBC(int,int,int,int);
    
        // methods
        virtual BoundaryCondition *PrintBC(ostream &);
        MatPtLoadBC *AddMPFluxBC(double);
	
		// accessors
		virtual void SetBCValue(double);

};

extern MatPtHeatFluxBC *firstHeatFluxPt;

#endif
