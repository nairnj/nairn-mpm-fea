/********************************************************************************
    NodalVelBC.hpp
    NairnMPM
    
    Created by John Nairn on Thu Apr 1 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
	
	Dependencies
		BoundaryCondition.hpp
********************************************************************************/

#ifndef _NODALVELBC_

#define _NODALVELBC_

#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Nodes/NodalPoint.hpp"

class NodalPoint;

class NodalVelBC : public BoundaryCondition
{
    public:
        int dir;
		Vector *pk;
 		double skewAngle;
        
        // constructors and destructors
        NodalVelBC(int,int,int,double,double);
		virtual ~NodalVelBC();
		virtual BoundaryCondition *UnsetDirection(void);
		virtual BoundaryCondition *SetRigidProperties(int,int,int,double);
        
        // virtual methods
        virtual BoundaryCondition *PrintBC(ostream &);
				
		// specific methods
        NodalVelBC *CopyNodalVelocities(NodalPoint *);
        NodalVelBC *PasteNodalVelocities(NodalPoint *);
		
		// class methods
		static void GridMomentumConditions(int);
		static void ConsistentGridForces(void);
		
		// accessors
		void SetSkewAngle(double);
};

// variables (changed in MPM time step)
extern NodalVelBC *firstVelocityBC;
extern NodalVelBC *lastVelocityBC;
extern NodalVelBC *firstRigidVelocityBC;
extern NodalVelBC *reuseRigidVelocityBC;

#endif

