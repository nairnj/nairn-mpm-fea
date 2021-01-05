/********************************************************************************
    LinearTraction.cpp
    nairn-mpm-fea
    
    Created by John Nairn on 4/1/08.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/LinearTraction.hpp"
#include "Cracks/CrackSegment.hpp"
#include "System/UnitsController.hpp"

#pragma mark LinearTraction::Constructors and Destructors

// Constructor 
LinearTraction::LinearTraction(char *matName,int matID) : CohesiveZone(matName,matID)
{
	kI1=kII1=0.;
}

#pragma mark LinearTraction::Initialization

// Check for required  slopes
const char *LinearTraction::VerifyAndLoadProperties(int np)
{
	// must always non-negative k1
	if(kI1<0. || kII1<0.)
		return "Linear traction law requires kIe>0 and kIIe>0";
		
	return TractionLaw::VerifyAndLoadProperties(np);
}

// print to output window
void LinearTraction::PrintMechanicalProperties(void) const
{
	PrintProperty("kI",kI1*UnitsController::Scaling(1.e-6),UnitsController::Label(TRACTIONSLOPE_UNITS));
	PrintProperty("kII",kII1*UnitsController::Scaling(1.e-6),UnitsController::Label(TRACTIONSLOPE_UNITS));
    cout << endl;
}

// Doesn't need parent class history
char *LinearTraction::InitHistoryData(char *pchr) { return NULL; }

#pragma mark LinearTraction::Traction Law

// Traction law - assume trianglar shape with unloading from down slope back to the origin
void LinearTraction::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,Vector *n,Vector *t,double area)
{
	double Tn=0.,Tt=0.;
	
	// normal force (only if open)
	if(nCod>0.)
		Tn = kI1*nCod;
	
	// shear (always)
	Tt = kII1*tCod;
	
	// force is traction times area projected onto plane of unit vectors (units F)
	// tract = -area*(Tn*n + Tt*t)
	// In 2D, if t=(dx,dy), then n=(-dy,dx)
	cs->tract.x = -area*(Tn*n->x + Tt*t->x);
	cs->tract.y = -area*(Tn*n->y + Tt*t->y);
	cs->tract.z = -area*(Tn*n->z + Tt*t->z);
}

// Return current traction law strain energy (Int T.du)
//	This energy is needed for J integral (and only used in J Integral)
// units of F/L
double LinearTraction::CrackWorkEnergy(CrackSegment *cs,double nCod,double tCod)
{
	// normal only if opend
	double workEnergy = nCod>0. ? 0.5*kI1*nCod*nCod : 0.;

	// add sheear energy
	workEnergy += 0.5*kII1*tCod*tCod;
	
	return workEnergy;
}

#pragma mark LinearTraction::Accessors

// return material type
const char *LinearTraction::MaterialType(void) const { return "Linear Elastic Traction"; }


