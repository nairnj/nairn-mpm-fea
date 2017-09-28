/********************************************************************************
	AdhesionFriction.cpp
	nairn-mpm-fea

	Friction sliding with adhesion and velocity-dependent coefficient of
	friction (linear dependence)

	Created by John Nairn, Nov 16, 2015.
	Copyright (c) 2015 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/AdhesionFriction.hpp"
#include "System/UnitsController.hpp"
#include "Cracks/CrackSurfaceContact.hpp"

#pragma mark AdhesionFriction::Constructors and Destructors

// Constructors
AdhesionFriction::AdhesionFriction() {}

AdhesionFriction::AdhesionFriction(char *matName) : CoulombFriction(matName)
{
	Sa = 0.;
	Na = 0.;
	kmu = 0.;
	vhalf = 0.;
	smoothStaticToDynamic = false;
}

#pragma mark AdhesionFriction::Initialization

// Read material properties by name (in xName). Set input to variable type
// (DOUBLE_NUM or INT_NUM) and return pointer to the class variable
// (cast as a char *)
char *AdhesionFriction::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    if(strcmp(xName,"Sa")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&Sa,gScaling,1.e6);
	}
	
    else if(strcmp(xName,"Na")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&Na,gScaling,1.e6);
	}
	
    else if(strcmp(xName,"kmu")==0)
	{	//  Legacy sec/mm
		input=DOUBLE_NUM;
		return (char *)&kmu;
	}
	
    else if(strcmp(xName,"vhalf")==0)
	{	//  Legacy mm/se
		input=DOUBLE_NUM;
		smoothStaticToDynamic = true;
		return (char *)&vhalf;
	}
	
	// does not all any from MaterialBase
    return CoulombFriction::InputMaterialProperty(xName,input,gScaling);
}

// Verify input properties do calculations; if problem return string with an error message
// Don't pass on to material base
const char *AdhesionFriction::VerifyAndLoadProperties(int np)
{
	if(frictionCoeff<0.)
		return "The coefficient of friction must be nonnegative";
	if(Sa<0. || Na<0.)
		return "The adhesive strengths must be nonnegative";
	
	// must call super class
	const char *err = CoulombFriction::VerifyAndLoadProperties(np);
	
	// reset if was changed
	frictionStyle = FRICTIONAL;
	
	// skip smooth transition if not needed
	if(vhalf<=0.) smoothStaticToDynamic = false;
	if(frictionCoeffStatic<=frictionCoeff) smoothStaticToDynamic = false;
	if(!smoothStaticToDynamic) vhalf = 0.;
	
	return err;
}

// print contact law details to output window
void AdhesionFriction::PrintContactLaw(void) const
{
	char hline[200];
	
	sprintf(hline,"Contact by friction and adhesion with coefficient of friction: %.6f",frictionCoeff);
	cout << hline << endl;
	
	if(frictionCoeffStatic>0.)
	{	sprintf(hline,"                           and static coefficient of friction: %.6f",frictionCoeffStatic);
		cout << hline << endl;
	}
	
	const char *label = UnitsController::Label(PRESSURE_UNITS);
	PrintProperty("Sa",Sa*UnitsController::Scaling(1.e-6),label);
	PrintProperty("Na",Na*UnitsController::Scaling(1.e-6),label);
	cout << endl;
	
	sprintf(hline,"1/(%s)",UnitsController::Label(CUVELOCITY_UNITS));
	PrintProperty("kmu",kmu,hline);
	PrintProperty("vhalf",vhalf,UnitsController::Label(CUVELOCITY_UNITS));
	cout << endl;

}

#pragma mark AdhesionFriction:Step Methods

// Return Sslide Ac dt = f(N) Ac dt
// Input is N Ac dt (and is always positive when in contact)
// If needed in the friction law, Ac is the contact area
// The relative sliding speed after correcting the momentum will be (SStickAcDt-SslideAcDt)/mred
double AdhesionFriction::GetSslideAcDt(double NAcDt,double SStickAcDt,double mred,
									   double contactArea,bool &inContact,double deltime) const
{
	// check for adhesion not exceeded
	if(!inContact)
	{	// if not failed, then it sticks (needs contact area)
		// note: never gets here unless Sa and Na are > 0
		double shear = SStickAcDt/(Sa*contactArea*deltime);
		double normal = NAcDt/(Na*contactArea*deltime);
		double surface = shear*shear + normal*normal;
		if(surface<1.)
		{	// adhesion not exceeded, so make it stick
			inContact = true;
			return SStickAcDt+1.;
		}
		
		// adhesion failure: leave inContact as false and return with any number
		return 0.;
	}

	// if in contact, return frictional sliding
	// Sslide = Sa + mu N or Sslide Ac dt = mu N Ac dt + Sa Ac dt
	// and mu = mu0 + kmu*vrel
	double Sslide = Sa>0. ? Sa*contactArea*deltime : 0. ;
	if(frictionCoeffStatic>0.)
	{	double Stest = Sslide + frictionCoeffStatic*NAcDt;
		if(SStickAcDt<Stest) return Stest;			// it will stick
	}
	
	// add friction terms as needed needed
	if(smoothStaticToDynamic)
	{	// smooth static dynamicand optionally adhesion and velocity dependent dynamic
		double kprime = kmu/mred;
		double Aterm = 1.+kprime*NAcDt;
		double Bterm = frictionCoeff + kprime*SStickAcDt + Sslide/NAcDt;
		double kstar = vhalf*mred;
		double Dstar = SStickAcDt + kstar;
		double arg0 = Bterm*NAcDt/Aterm;
		double arg1 = 0.5*(Dstar + arg0);
		double arg2 = 0.5*(Dstar - arg0);
		Sslide = arg1 - sqrt(arg2*arg2 - kstar*NAcDt*(frictionCoeffStatic-frictionCoeff)/Aterm);
		
	}
	else if(kmu!=0.)
	{	// sharp static/dynamic but linear velocity dependence and maybe adhesion
		double kmred = kmu/mred;
		Sslide = (Sslide + (frictionCoeff+kmred*SStickAcDt)*NAcDt)/(1. + kmred*NAcDt);
	}
	else
	{	// sharp static/dynamic and possibly adhesion
		Sslide += frictionCoeff*NAcDt;
	}
		
	// return result
	return Sslide;
}

#pragma mark AdhesionFriction::Accessors

// return unique, short name for this material
const char *AdhesionFriction::MaterialType(void) const { return "Adhesion and Friction"; }

// Continue to frictional calculation, even if not in contact, but only if has some separation adhesion
bool AdhesionFriction::ContactIsDone(bool inContact) const
{	if(inContact || (Sa>0. && Na>0.) ) return false;
	return true;
}

// needs area
bool AdhesionFriction::ContactLawNeedsContactArea(void) const { return Sa>0.; }


