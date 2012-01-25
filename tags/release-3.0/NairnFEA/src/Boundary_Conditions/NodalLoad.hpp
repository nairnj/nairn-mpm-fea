/********************************************************************************
    NodalLoad.hpp
    NairnFEA
    
    Created by John Nairn on Wed Mar 11 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _NODALLOAD_

#define _NODALLOAD_

#include "Boundary_Conditions/FEABoundaryCondition.hpp"

class NodalLoad : public FEABoundaryCondition
{
    public:

        // constructors and destructors
        NodalLoad(int,int);
        
        // methods
        NodalLoad *PrintLoad(void);
        NodalLoad *Reaction(double *,int,int);
		NodalLoad *MapNodes(int *);
		bool SameDofSetting(NodalLoad *);
};

extern NodalLoad *firstLoadBC;

#endif


