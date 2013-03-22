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
        
        // constructors and destructors
        NodalVelBC(int,int,int,double,double,double,double);
		virtual ~NodalVelBC();
		virtual BoundaryCondition *UnsetDirection(void);
		virtual BoundaryCondition *SetRigidProperties(int,int,int,double);
        
        // virtual methods
        virtual BoundaryCondition *PrintBC(ostream &);
				
		// specific methods
        NodalVelBC *CopyNodalVelocities(NodalPoint *);
        NodalVelBC *PasteNodalVelocities(NodalPoint *);
		NodalVelBC *ZeroVelBC(double);
		NodalVelBC *AddVelBC(double);
		NodalVelBC *SetGhostVelBC(double);
		NodalVelBC *InitFtot(double);
		NodalVelBC *AddFtot(double);
        int ConvertToDirectionBits(int);
        int ConvertToInputDof(void);
	
		// class methods
		static void GridMomentumConditions(int);
		static void ConsistentGridForces(void);
	
	protected:
		double currentValue;
        double angle1,angle2;
		Vector norm;
};

// variables (changed in MPM time step)
extern NodalVelBC *firstVelocityBC;
extern NodalVelBC *lastVelocityBC;
extern NodalVelBC *firstRigidVelocityBC;
extern NodalVelBC *reuseRigidVelocityBC;

#endif

