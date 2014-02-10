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
}

#pragma mark PressureLaw::Initialization

// no properties to read
char *PressureLaw::InputMat(char *xName,int &input)
{
    if(strcmp(xName,"stress")==0)
	{	input=DOUBLE_NUM;
		return((char *)&stress1);
	}
	
	else if(strcmp(xName,"function")==0)
	{	input=STRESS_FUNCTION_BLOCK;
		return((char *)this);
	}
	
    return TractionLaw::InputMat(xName,input);
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
	// Multiply by 1e6 to get N/mm/mm^2 (kg-m/sec^2/mm/mm^2) to g-mm/sec^2 / mm / mm^2
	sc = stress1*1.e6;
	
	// go to parent
	return TractionLaw::VerifyAndLoadProperties(np);
}

// print to output window
void PressureLaw::PrintMechanicalProperties(void) const
{
	if(function==NULL)
	{	PrintProperty("Stress",stress1,"MPa");
		cout <<  endl;
	}
	else
	{	char *expr=function->Expr('#');
		cout << "Stress = " << expr << endl;
	}
}

#pragma mark PressureLaw::Traction Law

// Traction law - constant pressure on crack surface
void PressureLaw::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,double dx,double dy,double area)
{
	double Tn;
	
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
	// physcial model is not as damage and therefor no unloading energy
	if(!fullEnergy) return 0.;
	
	// area at currnt time is stress*normal cod
	double tEnergy=0.;
	
	// constant pressure
	double tstress;
	if(function==NULL)
		tstress = stress1;
	else
	{	varTime = 1000.*mtime;
		tstress = function->Val();
	}
	
	// normal energy only if opened
	if(nCod>0.)
	{	tEnergy = tstress*nCod;		// N/mm
	}
	
	return tEnergy;
}

#pragma mark CubicTraction::Accessors

// Return the material tag
int PressureLaw::MaterialTag(void) const { return PRESSURELAWMATERIAL; }

// return material type
const char *PressureLaw::MaterialType(void) const { return "Pressure Traction Law"; }




