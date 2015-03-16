/********************************************************************************
	PressureLaw.cpp
	nairn-mpm-fea

	Created by John Nairn on 9/19/2013.
	Copyright (c) 2013 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/PressureLaw.hpp"
#include "Cracks/CrackSegment.hpp"
#include "Read_XML/mathexpr.hpp"

extern double mtime;

// global expression variables
double PressureLaw::varTime=0.;
PRVar plTimeArray[1] = { NULL };

#pragma mark PressureLaw::Constructors and Destructors

// Constructors with arguments
PressureLaw::PressureLaw(char *matName) : TractionLaw(matName)
{
	function = NULL;
	minCOD = -1.;
}

#pragma mark PressureLaw::Initialization

// no properties to read
char *PressureLaw::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    if(strcmp(xName,"stress")==0)
	{	input=DOUBLE_NUM;
		gScaling = 1.e6;			// Convert MPa to Pa
		return (char *)&stress1;
	}
	
	else if(strcmp(xName,"function")==0)
	{	input=STRESS_FUNCTION_BLOCK;
		return (char *)this;
	}
	
	else if(strcmp(xName,"minCOD")==0)
	{	input=DOUBLE_NUM;
		return (char *)&minCOD;
	}
	
    return TractionLaw::InputMaterialProperty(xName,input,gScaling);
}

// setting function if needed
void PressureLaw::SetStressFunction(char *bcFunction)
{	
	// NULL or empty is an error
	if(bcFunction==NULL)
		ThrowSAXException("Stress setting function of time is missing");
	if(strlen(bcFunction)==0)
		ThrowSAXException("Stress setting function of time is missing");
	
	// create time variable if needed
	if(plTimeArray[0]==NULL)
	{	plTimeArray[0]=new RVar("t",&varTime);
	}
	
	if(function!=NULL)
		ThrowSAXException("Duplicate stress setting function");
	function=new ROperation(bcFunction,1,plTimeArray);
	if(function->HasError())
		ThrowSAXException("Stress setting function is not valid");
}

// calculate properties used in analyses - here constant pressure law
const char *PressureLaw::VerifyAndLoadProperties(int np)
{
	// go to parent
	return TractionLaw::VerifyAndLoadProperties(np);
}

// print to output window
void PressureLaw::PrintMechanicalProperties(void) const
{
	if(function==NULL)
	{	PrintProperty("Stress",stress1*1.e-6,"MPa");
		cout <<  endl;
	}
	else
	{	char *expr=function->Expr('#');
		cout << "Stress = " << expr << endl;
	}
	if(minCOD>=0.)
	{	PrintProperty("Min COD",minCOD,"mm");
		cout << endl;
	}
}

#pragma mark PressureLaw::Traction Law

// Traction law - constant pressure on crack surface
void PressureLaw::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,double dx,double dy,double area)
{
	double Tn;
	
	// no pressure if less then a specific critical COD
	if(minCOD>=0. && nCod<=minCOD)
	{	cs->tract.x = 0.;
		cs->tract.y = 0.;
		return;
	}
	
	// constant pressure
	if(function==NULL)
		Tn = sc;
	else
	{	varTime = 1000.*mtime;
		Tn = function->Val()*1.e6;
	}
	
	// force is traction times area projected onto x-y plane
	cs->tract.x = area*Tn*dy;
	cs->tract.y = -area*Tn*dx;
}

// return total energy (which is needed for path independent J) under traction law curve
//		when fullEnergy is true
// return released enegery = total energy - recoverable energy (due to elastic unloading)
//		when fullEnergy is false
// units of N/mm
double PressureLaw::CrackTractionEnergy(CrackSegment *cs,double nCod,double tCod,bool fullEnergy)
{
	// physcial model is not as damage and therefore no unloading energy
	// also zero energy if close (nCod<=0)
	if(!fullEnergy || nCod<=0.) return 0.;
	
	// no pressure if less then a specific critical COD
	if(minCOD>=0. && nCod<=minCOD) return 0.;
	
	// constant pressure
	double tstress;
	if(function==NULL)
		tstress = stress1;
	else
	{	varTime = 1000.*mtime;
		tstress = function->Val();
	}
	
	return tstress*nCod;		// N/mm
}

#pragma mark CubicTraction::Accessors

// Return the material tag
int PressureLaw::MaterialTag(void) const { return PRESSURELAWMATERIAL; }

// return material type
const char *PressureLaw::MaterialType(void) const { return "Pressure Traction Law"; }




