/********************************************************************************
    NodalVelBC.hpp
    nairn-mpm-fea
    
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
		static bool holdAllVelocityBCs;
        
        // constructors and destructors
        NodalVelBC(int,int,int,double,double,double,double);
        virtual int GetSetDirection(void) const;
		virtual BoundaryCondition *SetRigidProperties(int,int,int,double);
        
        // virtual methods
        virtual BoundaryCondition *PrintBC(ostream &);
				
		// boundary conditions
		virtual NodalVelBC *GetCurrentBCValue(double);
		virtual NodalVelBC *ZeroVelocityBC(double,int);
		virtual NodalVelBC *AddVelocityBC(double,int);
	
		// other methods
		NodalVelBC *SetMirroredVelBC(double);
		NodalVelBC *AddReactionForce(Vector *,int);
        int ConvertToDirectionBits(int);
        void SetNormalVector(void);
		Vector *GetNormalVector(void);
        int ConvertToInputDof(void);
		void SetReflectedNode(int,double);
		void SetMirrorSpacing(int);
	
		// class methods
		static void GridVelocityBCValues(void);
		static void GridVelocityConditions(int);
		static void VelocityBCLoop(int);
		static Vector TotalReactionForce(int);
	
	protected:
		double currentValue;
        double angle1,angle2;
		Vector norm;
		Vector freaction;
        int mirrorSpacing;
		int reflectedNode;
		double reflectRatio;
};

// variables (changed in MPM time step)
extern NodalVelBC *firstVelocityBC;
extern NodalVelBC *lastVelocityBC;
extern NodalVelBC *firstRigidVelocityBC;
extern NodalVelBC *reuseRigidVelocityBC;

#endif

