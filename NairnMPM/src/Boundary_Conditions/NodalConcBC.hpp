/********************************************************************************
    NodalConcBC.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Apr 1 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
	
	Dependencies
		BoundaryCondition.hpp
********************************************************************************/

#ifndef _NODALCONCBC_

#define _NODALCONCBC_

#include "Boundary_Conditions/NodalValueBC.hpp"
#include "Materials/MaterialBase.hpp"

class NodalPoint;

class NodalConcBC: public NodalValueBC
{
    public:
		int phaseStyle;
    
        // constructors and destructors
        NodalConcBC(int,int,double,double,int);
	
		// methods
		virtual BoundaryCondition *PrintBC(ostream &);
	
        // Accessors
		virtual void SetBCValue(double bcvalue);
        virtual int GetSetDirection(void) const;
};

// variables (changed in MPM time step)
extern NodalConcBC *firstConcBC;
extern NodalConcBC *lastConcBC;
extern NodalConcBC *firstRigidConcBC;
extern NodalConcBC *reuseRigidConcBC;
extern NodalConcBC *firstDiffBC[NUM_DUFFUSION_OPTIONS];

#endif

