/********************************************************************************
	SmoothStep3.cpp
	nairn-mpm-fea
 
	Softening law is linear f(delta) = 1 - S1[delta/deltamax]
 
    where S1[x] = -2 x^3 + 3 x^2 is the cubic smooth step function. This
    function has properties f(0)=1, f(deltamax) = 0, f'(0)=0, and
    f'(deltamax) = 0
 
 		s = gscaling = Ac/(Vp rho sigma(sp))
 
	Created by John Nairn, Dec 25, 2016.
	Copyright (c) 2016 John A. Nairn, All rights reserved.
 *******************************************************************************/

#include "stdafx.h"
#include "Materials/SmoothStep3.hpp"
#include "Materials/MaterialBase.hpp"
#include "System/UnitsController.hpp"

#pragma mark SmoothStep3::Constructors and Destructors

// Constructors
SmoothStep3::SmoothStep3() : SofteningLaw()
{
	// default to zero initial stiffness
	k = 0.;
	stepOnly = true;
}

#pragma mark SmoothStep3::Initialize

// Read softening law properties
char *SmoothStep3::InputSofteningProperty(char *xName,int &input,double &gScaling)
{
	// Fracture toughness
	if(strcmp(xName,"k")==0)
	{	input=DOUBLE_NUM;
		return (char *)&k;
	}
	
	return SofteningLaw::InputSofteningProperty(xName,input,gScaling);
}

// print just initiation properties to output window
void SmoothStep3::PrintSofteningProperties(double sigmainit)
{
	// print initial stiffness info (if was entered)
	double pkRatio = 1.;
	if(k>0)
	{	k2 = 1.+0.5*k;
		double k3 = 1.+k/3.;
		k6 = 1.+k/6.;
		pkRatio = (k2/k3)*(k2/k3)/k3;
		kterm = (4.*k2-1.)/3.;
		stepOnly = false;
	}
	else
	{	k = 0.;
		k2 = 1.;
		kterm = 1.;
		k6 = 1.;
	}
	
	// prints namem Gc, ucrit, and ecrti
	SofteningLaw::PrintSofteningProperties(sigmainit);
	
	// Initial stiffness
	MaterialBase::PrintProperty("kinit",k,"");
	
	// peak stress
	double sigmapeak = sigmainit/pkRatio;
	MaterialBase::PrintProperty("s(max)",sigmapeak*UnitsController::Scaling(1.e-6),"");
	cout << endl;
}

#pragma mark SmoothStep3::Required Methods

// initiation law name
const char *SmoothStep3::GetSofteningLawName(void) const { return "Cubic step function softening"; }

// f(delta) relative to initiation stress
// gScaling must include relative toughness
double SmoothStep3::GetFFxn(double delta,double gScaling) const
{	double deltaMax = 2.*Gc*gScaling/k6;
    double x = delta/deltaMax;
	double arg = 1.-x;
    return (1.+2.*k2*x)*arg*arg;
}

// gScaling must include relative toughness
double SmoothStep3::GetDeltaMax(double gScaling) const { return 2.*Gc*gScaling/k6; }

// f'(delta)
// gScaling must include relative toughness
double SmoothStep3::GetFpFxn(double delta,double gScaling) const
{	double deltaMax = 2.*Gc*gScaling/k6;
    double x = delta/deltaMax;
    return -(6.*k2*x-k)*(1.-x)/deltaMax;
}

// Get energy released (per unit volume per unit stress or Gbar/sigma) up to delta.
// gScaling must include relative toughness
double SmoothStep3::GetGToDelta(double delta,double gScaling) const
{	double deltaMax = 2.*Gc*gScaling/k6;
    double x = delta/deltaMax;
	double x2 = x*x;
    return 0.5*delta*(1. + x2*(kterm - k2*x));
}

// Get G/Gc up to delta or Gbar(delta)/Gbar(infinity)
// gScaling must include relative toughness
double SmoothStep3::GetGoverGc(double delta,double gScaling) const
{	double deltaMax = 2.*Gc*gScaling/k6;
    double x = delta/deltaMax;
	double x2 = x*x;
    return x*(1. + x2*(kterm - k2*x))/k6;
}

// dimensionless stability factor = maxSlope * sGc
double SmoothStep3::GetEtaStability(void) const
{ 	if(stepOnly) return 4./3.;
	double otk2 = (1.+2*k2);
	return 12.*k2/(otk2*otk2*k6);
}

// Find Phi(delta) function = f(delta)-delta*f'(delta))
// gscaling must include relative values
double SmoothStep3::GetPhiFxn(double delta,double gScaling) const
{	double deltaMax = 2.*Gc*gScaling/k6;
	double x = delta/deltaMax;
	double x2 = x*x;
	return 1. + x2*(4*k2*(1.-x)-1.);
}

// Find Rd(delta) function = ei(f(delta)-delta*f'(delta))/(delta+ei f(delta)_^2
// gscaling and e0 must include relative values
double SmoothStep3::GetRdFxn(double delta,double gScaling,double e0) const
{	double deltaMax = 2.*Gc*gScaling/k6;
	double x = delta/deltaMax;
	double x2 = x*x;
	double arg = 1.-x;
	double en = delta + e0*((1.+2.*k2*x)*arg*arg);
	return e0*(1. + x2*(4*k2*(1.-x)-1.))/(en*en);
}

#pragma mark SmoothStep3::Optional Methods

#ifdef SS3_ANALYTICAL
// Solve for increment in crack opening strain by solving de = ddelta + e0*(f(delta+ddelta)-f(delta))
// Return ddelta if not failed or -1 if failed
// gScaling must include relative toughness, e0 must include relative strength
double SmoothStep3::GetDDelta(double de,double e0,double delta,double d,double gScaling) const
{	double deltaMax = 2.*Gc*gScaling/k6;
    double x = delta/deltaMax;
	double arg,arg2;
	if(stepOnly)
	{	arg = (1.-2.*x)*deltaMax - 2.*de - (1.+4.*x*x*(x-1.5))*e0;
		arg2 = 2.*deltaMax - 3.*e0;
	}
	else
	{	double k22 = k2*k2;
		double k23 = k22*k2;
		arg = k2*(kterm-2.*k2*x)*deltaMax - 2.*k2*k2*de
				- e0*((-1.-6.*k2+42.*k22-8.*k23)/27. +4.*k22*(k2-1.)*x + 4.*k22*x*x*(k2*x-0.5*(4.*k2-1.)));
		double otk2 = 1.+2.*k2;
		arg2 = 2.*k2*deltaMax - otk2*otk2*e0/3.;
	}
	double denom = pow(9.*e0*e0*arg + e0*sqrt(3.*e0*(arg2*arg2*arg2 + 27.*e0*arg*arg)),1./3.);
    double ddelta = 0.5*deltaMax*( kterm/k2- 2.*x + arg2/(CUBEROOT3*k2*denom) - denom/(CUBEROOT32*k2*e0) );
    
    // has it failed
    if(delta+ddelta>deltaMax) ddelta = -1.;
    return ddelta;
}
#endif


