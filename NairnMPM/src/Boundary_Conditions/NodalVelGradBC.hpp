/********************************************************************************
 NodalVelGradBC.hpp
 nairn-mpm-fea
 
 Created by John Nairn on Feb 27, 2019.
 Copyright (c) 2019 John A. Nairn, All rights reserved.
 
 Dependencies
 	NodalVelBC.hpp, BoundaryCondition.hpp, NodalPoint.hpp
 ********************************************************************************/

#ifndef _NODALVELGRADBC_

#define _NODALVELGRADBC_

#include "Boundary_Conditions/NodalVelBC.hpp"

class NodalVelGradBC : public NodalVelBC
{
	public:
	
		// constructors and destructors
		NodalVelGradBC(int,int,int,double,double);
		virtual ~NodalVelGradBC();
	
		// methods
		virtual BoundaryCondition *PrintBC(ostream &);
		virtual NodalVelGradBC *GetCurrentBCValue(double);
		virtual NodalVelGradBC *ZeroVelocityBC(double,int);
		virtual NodalVelGradBC *AddVelocityBC(double,int);
	
		// accessors
		int ConvertToDirectionBits(int);
		void SetGradFunction(char *,bool);
		void BCGradValue(double,int,int);
		void BCPosition(double,int);

	protected:
		int side;
		double depth;
		Expression *gradFunction;
		Expression *dispFunction;
		double gradValue;				// scalars only controlled axis
		double position;
		int firstNode,nodeStep;
		bool extractGrad;			// true when grad function is mirror

};

#endif
