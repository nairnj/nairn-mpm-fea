/********************************************************************************
    NodalConcBC.hpp
    NairnMPM
    
    Created by John Nairn on Thu Apr 1 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
	
	Dependencies
		BoundaryCondition.hpp
********************************************************************************/

#ifndef _NODALCONCBC_

#define _NODALCONCBC_

#include "Boundary_Conditions/BoundaryCondition.hpp"

class NodalPoint;

class NodalConcBC: public BoundaryCondition
{
    public:
		double concentrationNoBC;
        
        // constructors and destructors
        NodalConcBC(int,int,double,double);
		virtual BoundaryCondition *UnsetDirection(void);
		virtual BoundaryCondition *SetRigidProperties(int,int,int,double);
        
        // methods
        BoundaryCondition *PrintBC(ostream &);
		NodalConcBC *CopyNodalConcentration(NodalPoint *);
		NodalConcBC *PasteNodalConcentration(NodalPoint *);

};

// variables (changed in MPM time step)
extern NodalConcBC *firstConcBC;
extern NodalConcBC *lastConcBC;
extern NodalConcBC *firstRigidConcBC;
extern NodalConcBC *reuseRigidConcBC;

#endif

