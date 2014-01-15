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

#include "Boundary_Conditions/BoundaryCondition.hpp"

class NodalPoint;

class NodalTempBC: public BoundaryCondition
{
    public:
		double temperatureNoBC;
        
        // constructors and destructors
        NodalTempBC(int,int,double,double);
		virtual BoundaryCondition *UnsetDirection(void);
		virtual BoundaryCondition *SetRigidProperties(int,int,int,double);
        
        // methods
        BoundaryCondition *PrintBC(ostream &);
		NodalTempBC *CopyNodalTemperature(NodalPoint *);
		NodalTempBC *PasteNodalTemperature(NodalPoint *);

};

// variables (changed in MPM time step)
extern NodalTempBC *firstTempBC;
extern NodalTempBC *lastTempBC;
extern NodalTempBC *firstRigidTempBC;
extern NodalTempBC *reuseRigidTempBC;

#endif

