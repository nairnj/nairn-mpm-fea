/********************************************************************************
    MatPtLoadBC.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Tues Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Dependencies
		BoundaryCondition.hpp
********************************************************************************/

#ifndef _MATPTLOADBC_

#define _MATPTLOADBC_

#include "Boundary_Conditions/BoundaryCondition.hpp"

class MatPtLoadBC: public BoundaryCondition
{
    public:
        int ptNum;
        int direction;
        
        // constructors and destructors
        MatPtLoadBC(int,int,int);
		MatPtLoadBC *ReorderPtNum(int,int);
    
        // virtual methods
        virtual BoundaryCondition *PrintBC(ostream &);
	
		// specific methods (not overridden)
        MatPtLoadBC *ZeroMPLoad(void);
        MatPtLoadBC *AddMPLoad(double);
        MatPtLoadBC *ReverseLinearLoad(double,double *,bool);
        MatPtLoadBC *MakeConstantLoad(double);
	
		// accessors
		virtual void GetPosition(unordered_map<string,double>);
		virtual void SetBCValue(double);
    
		// overridden by flux sub classes
		virtual MatPtLoadBC *AddMPFluxBC(double);
	
		// class methods
		static void SetParticleFext(double);
    
    private:
        double holdValue;
        bool holding;
};

extern MatPtLoadBC *firstLoadedPt;

#endif
