/********************************************************************************
    InitialCondition.cpp
    nairn-mpm-fea

    Created by John Nairn on July 5, 2017
    Copyright (c) 2017 John A. Nairn, All rights reserved.
 ********************************************************************************/

#include "stdafx.h"
#include "Boundary_Conditions/InitialCondition.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Materials/MaterialBase.hpp"

// global
InitialCondition *firstDamagedPt=NULL;

#pragma mark InitialCondition: Constructors and Destructors

InitialCondition::InitialCondition(int cType,int num) : MatPtLoadBC(num,1,1)
{
    icType = cType;
    dnorm = MakeVector(1.,0.,0.);
    dvals = MakeVector(0.,0.,0.);
}

#pragma mark InitialCondition: Methods

InitialCondition *InitialCondition::AssignInitialConditions(bool is3D)
{
    int matid=mpm[ptNum-1]->MatID();
    if(matid>=nmat)
        throw CommonException("Damaged material point with an undefined material type","InitialCondition::AssignInitialConditions");
    
    // Membrane particle need membrane material and solid particle need solid material
    if(theMaterials[matid]->MaterialStyle()!=SOLID_MAT)
        throw CommonException("Damaged material point with non-solid material","InitialCondition::AssignInitialConditions");
    
    // initial conditions on this particle
    theMaterials[matid]->SetInitialConditions(this,mpm[ptNum-1],is3D);
   
    return (InitialCondition *)GetNextObject();
}

#pragma mark InitialCondition: Accessors

// set vector value
void InitialCondition::SetInitialDamage(Vector *nvec,Vector *dvec,double mode)
{
    // normal (be sure it is normalized)
    dnorm = *nvec;
    double scale = sqrt(dnorm.x*dnorm.x + dnorm.y*dnorm.y + dnorm.z*dnorm.z);
    dnorm.x /= scale;
    dnorm.y /= scale;
    dnorm.z /= scale;
    
    // three damage variables (expand if ever needed)
    dvals = *dvec;
	
	// damage mode for history [1]
	dmode = mode;
}

// get damage vectors for normal and for parameters
Vector InitialCondition::GetDamageNormal(void) { return dnorm; };
Vector InitialCondition::GetDamageParams(double &mode)
{	mode = dmode;
	return dvals;
};


