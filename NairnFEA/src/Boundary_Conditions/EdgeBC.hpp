/********************************************************************************
    EdgeBC.hpp
    NairnFEA
    
    Created by John Nairn on Mon Mar 15 2004.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _EDGEBC_

#define _EDGEBC_

#include "Boundary_Conditions/FEABoundaryCondition.hpp"

class EdgeBC : public FEABoundaryCondition
{
    public:
       // constructors and destructors
        EdgeBC(int,int,int);
        
        // methods
		EdgeBC *PrintEdgeLoad(void);
		void GetConsistentLoads(double *,int);
		void SetStress(double *,int);
		bool SameDofSetting(EdgeBC *);
		int ElementIndex(void);
	
	private:
		int element;
		int face,facenodes;
		double stress[3];
};

extern EdgeBC *firstEdgeBC;

#endif
