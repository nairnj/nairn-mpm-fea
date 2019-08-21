/********************************************************************************
    NodalTempBC.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Oct 18 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
	
	Dependencies
		BoundaryCondition.hpp
********************************************************************************/

#ifndef _NODALTEMPBC_

#define _NODALTEMPBC_

#include "Boundary_Conditions/NodalValueBC.hpp"

class NodalPoint;

class NodalTempBC: public NodalValueBC
{
    public:
		double *temperatureNoBC;
        
        // constructors and destructors
        NodalTempBC(int,int,double,double);
    
        // methods
 		NodalTempBC *CopyNodalValue(NodalPoint *);
		NodalTempBC *PasteNodalValue(NodalPoint *);
	
        // reaction heat calculation
		NodalTempBC *AddHeatReaction(double *,int);
	
		// class methods
		static double TotalHeatReaction(int);
    
        // accessors
        virtual int GetSetDirection(void) const;
        virtual TransportField *GetTransportFieldPtr(NodalPoint *) const;
};

// variables (changed in MPM time step)
extern NodalTempBC *firstTempBC;
extern NodalTempBC *lastTempBC;
extern NodalTempBC *firstRigidTempBC;
extern NodalTempBC *reuseRigidTempBC;

#endif

