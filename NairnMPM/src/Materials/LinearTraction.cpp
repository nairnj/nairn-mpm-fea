/********************************************************************************
    LinearTraction.cpp
    nairn-mpm-fea
    
    Created by John Nairn on 4/1/08.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/LinearTraction.hpp"
#include "Cracks/CrackSegment.hpp"
#include "System/UnitsController.hpp"

#pragma mark LinearTraction::Constructors and Destructors

// Constructors with arguments 
LinearTraction::LinearTraction(char *matName) : CohesiveZone(matName)
{
	kI1=kII1=0.;
}

#pragma mark LinearTraction::Initialization

// Check for required  slopes
const char *LinearTraction::VerifyAndLoadProperties(int np)
{
	// must always non-negative k1
	if(kI1<0. || kII1<0.)
		return "Linear traction law requires non-negative kIe and kIIe";
		
	return TractionLaw::VerifyAndLoadProperties(np);
}

// print to output window
void LinearTraction::PrintMechanicalProperties(void) const
{
	PrintProperty("kI",kI1*UnitsController::Scaling(1.e-6),UnitsController::Label(TRACTIONSLOPE_UNITS));
	PrintProperty("kII",kII1*UnitsController::Scaling(1.e-6),UnitsController::Label(TRACTIONSLOPE_UNITS));
    cout << endl;
}

#pragma mark LinearTraction::Traction Law

// Traction law - assume trianglar shape with unloading from down slope back to the origin
void LinearTraction::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,double dx,double dy,double area)
{
	double Tn=0.,Tt=0.;
	
	// normal force (only if open)
	if(nCod>0.)
		Tn = kI1*nCod;
	
	// shear (always)
	Tt = kII1*tCod;
	
	// force is traction time area projected onto x-y plane
	cs->tract.x = area*(Tn*dy - Tt*dx);
	cs->tract.y = area*(-Tn*dx - Tt*dy);
}

// return total energy (which is needed for path independent J) under traction law curve
//		when fullEnergy is true
// return released enegery = total energy - recoverable energy (due to elastic unloading)
//		when fullEnergy is false
// units of N/mm
double LinearTraction::CrackTractionEnergy(CrackSegment *cs,double nCod,double tCod,bool fullEnergy)
{
	// all the energy is recoverable since it is elastic
	if(!fullEnergy) return 0.;
	
	double tEnergy = 0.;
	
	// get entire area under the curve
	
	// normal energy only if opened
	if(nCod>0.)
	{	double Tn = kI1*nCod;
		tEnergy = 0.5*Tn*nCod;
	}
	
	double Tt = kII1*tCod;
	tEnergy += 0.5*Tt*tCod;
	
	return tEnergy;
}

#pragma mark LinearTraction::Accessors

// return material type
const char *LinearTraction::MaterialType(void) const { return "Linear Elastic Traction"; }

// Return the material tag
int LinearTraction::MaterialTag(void) const { return LINEARTRACTIONMATERIAL; }


