/********************************************************************************
	PressureLaw.cpp
	nairn-mpm-fea

	Created by John Nairn on 9/19/2013.
	Copyright (c) 2013 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/PressureLaw.hpp"
#include "Cracks/CrackSegment.hpp"
#include "System/UnitsController.hpp"
#include "Read_XML/Expression.hpp"

extern double mtime;

#pragma mark PressureLaw::Constructors and Destructors

// Constructor
PressureLaw::PressureLaw(char *matName,int matID) : TractionLaw(matName,matID)
{
	function = NULL;
	minCOD = -1.;
}

#pragma mark PressureLaw::Initialization

// no properties to read
char *PressureLaw::InputTractionLawProperty(char *xName,int &input,double &gScaling)
{
    if(strcmp(xName,"stress")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&stress1,gScaling,1.e6);
	}
	
	else if(strcmp(xName,"function")==0)
	{	input=STRESS_FUNCTION_BLOCK;
		return (char *)this;
	}
	
	else if(strcmp(xName,"minCOD")==0)
	{	input=DOUBLE_NUM;
		return (char *)&minCOD;
	}
	
    return TractionLaw::InputTractionLawProperty(xName,input,gScaling);
}

// setting function if needed
// throws std::bad_alloc, SAXException()
void PressureLaw::SetStressFunction(char *bcFunction)
{	
	// NULL or empty is an error
	if(bcFunction==NULL)
	{	ThrowSAXException("Stress setting function of time is missing");
		return;
	}
	if(strlen(bcFunction)==0)
	{	ThrowSAXException("Stress setting function of time is missing");
		return;
	}

	// repeat is an error
	if(function!=NULL)
		ThrowSAXException("Duplicate stress setting function");
	
	// create function
	function = Expression::CreateExpression(bcFunction,"Stress setting function is not valid");
}

// print to output window
void PressureLaw::PrintMechanicalProperties(void) const
{
	if(function==NULL)
	{	PrintProperty("Stress",stress1*UnitsController::Scaling(1.e-6),UnitsController::Label(PRESSURE_UNITS));
		cout <<  endl;
	}
	else
	{	cout << "Stress = " << function->GetString() << endl;
	}
	if(minCOD>=0.)
	{	PrintProperty("Min COD",minCOD,UnitsController::Label(CULENGTH_UNITS));
		cout << endl;
	}
}

// evaluate pressure at current time. If no function use input constant in stress 1
// (Don't call in parallel code due to function)
void PressureLaw::CalculateTimeFunction(void)
{
	if(function!=NULL)
	{	// in Legacy, convert MPa to Pa, in consistent units use the function
		stress1 = function->TValue(mtime*UnitsController::Scaling(1000.))*UnitsController::Scaling(1.e6);
	}
}

#pragma mark PressureLaw::Traction Law

// Traction law - constant pressure on crack surface
void PressureLaw::CrackTractionLaw(CrackSegment *cs, double nCod, double tCod, Vector *n, Vector *t, double area)
{
	// no pressure if less then a specific critical COD
	if(minCOD >= 0. && nCod <= minCOD)
	{
		cs->tract.x = 0.;
		cs->tract.y = 0.;
		cs->tract.z = 0.;
		return;
	}

	// force is traction times area projected onto plane of unit vectors (units F)
	// tract = -area*(Tn*n + Tt*t), in 2D, if t=(dx,dy), then n=(-dy,dx)
	cs->tract.x = -area*stress1*n->x;
	cs->tract.y = -area*stress1*n->y;
	cs->tract.z = -area*stress1*n->z;
}

// return total energy (which is needed for path independent J) under traction law curve
//		when fullEnergy is true
// return released enegery = total energy - recoverable energy (due to elastic unloading)
//		when fullEnergy is false
// units of N/mm
double PressureLaw::CrackTractionEnergy(CrackSegment *cs,double nCod,double tCod,bool fullEnergy)
{
	// physcial model is not as damage and therefore no unloading energy
	// also zero energy if closed (nCod<=0)
	if(!fullEnergy || nCod<=0.) return 0.;
	
	// no pressure if less then a specific critical COD
	if(minCOD>=0. && nCod<=minCOD) return 0.;
	
	// constant pressure
	return stress1*nCod;
}

#pragma mark CubicTraction::Accessors

// return material type
const char *PressureLaw::MaterialType(void) const { return "Pressure Traction Law"; }




