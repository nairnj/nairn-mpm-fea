/********************************************************************************
    InitialCondition.hpp
    nairn-mpm-fea

    Created by John Nairn on July 5, 2017
    Copyright (c) 2017 John A. Nairn, All rights reserved.
	
    Dependencies
        BoundaryCondition.hpp, MatPtLoad.hpp
********************************************************************************/

#ifndef _INITIALCONDITION_

#define _INITIALCONDITION_

#include "Boundary_Conditions/MatPtLoadBC.hpp"

// loading boundary conditions
enum {	UNKNOWN_CONDITION=0,INITIAL_DAMAGE };

class InitialCondition : public MatPtLoadBC
{
    public:
    
        // constructors
        InitialCondition(int,int);
    
        // methods
        virtual InitialCondition *AssignInitialConditions(bool);

        // accessors
        virtual void SetInitialDamage(Vector *,Vector *,double);
        Vector GetDamageNormal(void);
        Vector GetDamageParams(double &);
    
    protected:
        int icType;     // future use
        Vector dnorm,dvals;
		double dmode;
	
};

extern InitialCondition *firstDamagedPt;


#endif
