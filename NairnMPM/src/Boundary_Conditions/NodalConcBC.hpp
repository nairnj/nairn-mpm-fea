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

class NodalPoint;

class NodalConcBC: public NodalValueBC
{
    public:
    
        // constructors and destructors
        NodalConcBC(int,int,double,double);
	
		// methods
		virtual BoundaryCondition *PrintBC(ostream &);
	
        // Accessors
		virtual void SetBCValue(double bcvalue);
        virtual int GetSetDirection(void) const;
        virtual TransportField *GetTransportFieldPtr(NodalPoint *) const;
 
};

// variables (changed in MPM time step)
extern NodalConcBC *firstConcBC;
extern NodalConcBC *lastConcBC;
extern NodalConcBC *firstRigidConcBC;
extern NodalConcBC *reuseRigidConcBC;

#endif

