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

// Constructor
AdhesionFriction::AdhesionFriction(char *matName,int matID) : CoulombFriction(matName,matID)
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
char *AdhesionFriction::InputContactProperty(char *xName,int &input,double &gScaling)
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
    return CoulombFriction::InputContactProperty(xName,input,gScaling);
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
	
	// Do not allow stick or frictionless
	if(frictionStyle==STICK)
		return "AdhesiveFriction contact law does allow stick option";
	
	// reset if was changed to frictionless
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
	CoulombFriction::PrintContactLaw();
	
	char hline[200];
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
// Input is N Ac dt (and is always positive for compression)
// If needed in the friction law, Ac is the contact area
// The relative sliding speed after correcting the momentum will be (SStickAcDt-SslideAcDt)/mred
double AdhesionFriction::GetSslideAcDt(double NAcDt,double SStickAcDt,double mred,
									   double contactArea,bool &inContact,double deltime) const
{
	// check for adhesive failure when not in contact
	if(!inContact || NAcDt<0.)
	{	// If displacementOnly>0.1, then deln>0 and NAcDt may be + or -
		// otherwise either deln>0 or N<displacementOnly or both
		double surface = 2.;
		if(Sa>0. && Na>0.)
		{	double shear = SStickAcDt/(Sa*contactArea*deltime);
			// only look at tension (NAcDt<0)
			double normal = fmin(NAcDt,0.)/(Na*contactArea*deltime);
			surface = shear*shear + normal*normal;
		}
		
		// If surface<1 than should stick
		if(surface<1.)
		{	inContact = true;
			return SStickAcDt+1.;
		}
		
		// It has failed, exit if was determined not in contact
		if(!inContact) return 0;
		
		// Here means displacementOnly = 1 or < 0 with NAcDt>displacementOnly
		// Fall through to code, which accounts for possibly negative NAcDt
	}
	
	// Here inContact is true and NAcDt>do (which may be negative)
	double Sslide = Sa>0. ? Sa*contactArea*deltime : 0. ;
	
	// Return frictional sliding
	// Sslide = Sa + mu(eff) N or Sslide Ac dt = mu(eff) N Ac dt + Sa Ac dt
	// and mu(eff) = mu(trans) + kmu*vrel where mu(trans) is smooth static to dynamic transition
	
	// if static check for stick (here vrel=0)
	if(frictionCoeffStatic>0.)
	{	double Stest = Sslide + fmax(frictionCoeffStatic*NAcDt,0.);
		if(Stest>SStickAcDt) return Stest;			// it will stick
	}
	
	// or add velocity dependent friction terms as needed needed
	if(smoothStaticToDynamic && NAcDt>0.)
	{	// smooth static dynamic and optionally adhesion and velocity dependent dynamic
		// See notes on Contact Laws Paper and only solved for NAcDt>0.
		// If NAcDt<0 (i.e. displacmentOnly =1 or <0), ignores smooth transition
		double kstar = vhalf*mred;
		double Dstar = SStickAcDt + kstar;
		double arg0 = Sslide + frictionCoeff*NAcDt;
		double arg1 = 0.5*(Dstar + arg0);
		double arg2 = 0.5*(Dstar - arg0);
		double scale = 1.;
		double Aterm = 0.;
		if(kmu!=0.)
		{	double kmred = kmu/mred;
			arg1 += 0.5*kmred*NAcDt*(Dstar-kstar);
			scale = 1 + kmred*NAcDt;
			Aterm = kmred*(Dstar-kstar-Sslide-NAcDt*frictionCoeffStatic);
		}
		if(scale > 0.)
			Sslide = (arg1 - sqrt(arg2*arg2 + Aterm - kstar*NAcDt*(frictionCoeffStatic-frictionCoeff)))/scale;
		else if(kmu<0.)
		{	// stick if beyond the limit
			Sslide = SStickAcDt+1.;
		}
		else
		{	// slip if too low
			Sslide = 0.;
		}
	}
	else if(kmu!=0.)
	{	// sharp static/dynamic but linear velocity dependence and adhesion
		// Note: published paper as 1+kmed*dn = 1 -rmred*NAcDt, but that is a misprint. Sign correct here
		double kmred = kmu/mred;
		if(1.+kmred*NAcDt > 0.)
			Sslide = fmax((Sslide + (frictionCoeff+kmred*SStickAcDt)*NAcDt)/(1. + kmred*NAcDt),0.);
		else if(kmu<0.)
		{	// stick if beyond the limit
			Sslide = SStickAcDt+1.;
		}
		else
		{	// slip if too low
			Sslide = 0.;
		}
	}
	else
	{	// sharp static/dynamic with adhesion
		Sslide += fmax(frictionCoeff*NAcDt,0.);
	}
		
	// return result
	return Sslide;
}

#pragma mark AdhesionFriction::Accessors

// return unique, short name for this material
const char *AdhesionFriction::MaterialType(void) const { return "Adhesion and Friction"; }

// needs area
bool AdhesionFriction::ContactLawNeedsContactArea(void) const
{	if(Sa>0. && Na>0.) return true;
	return CoulombFriction::ContactLawNeedsContactArea();
}


