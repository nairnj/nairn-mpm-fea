/********************************************************************************
	ExponentialSoftening.cpp
	nairn-mpm-fea

	Softening law is linear f(delta) = exp(-k delta)

	where k depends on particle size and crack area by

 		1/k = (Ac/(Vp rho sigma(sp)) * Gc = s * Gc

	where s is scaling factor is input to functions that need it

 		s = gscaling = Ac/(Vp rho sigma(sp))

	Created by John Nairn, Dec 25, 2016.
	Copyright (c) 2016 John A. Nairn, All rights reserved.
*******************************************************************************/

#include "stdafx.h"
#include "Materials/ExponentialSoftening.hpp"

#pragma mark ExponentialSoftening::Constructors and Destructors

// none needed

#pragma mark ExponentialSoftening::Initialize

// none needed

#pragma mark ExponentialSoftening::Required Methods

// initiation law name
const char *ExponentialSoftening::GetSofteningLawName(void) const { return "Exponential softening"; }

// This law is expnential such that f(delta) = exp(-k*delta)
// gScaling must include relative toughness
double ExponentialSoftening::GetFFxn(double delta,double gScaling) const
{	double k = 1./(Gc*gScaling);
	return exp(-k*delta);
}

// This law has no maximum, but effective maximum is when numerical code
// is truncated at exp(-k dmax) = minFdelta or dmax = - ln(minFdelta)/k
// gScaling must include relative toughness
double ExponentialSoftening::GetDeltaMax(double gScaling) const
{	return unscaledDeltaMax*gScaling;
}

// This law is exponential such that such that f'(delta) = -k*exp(-k*delta)
// gScaling must include relative toughness
double ExponentialSoftening::GetFpFxn(double delta,double gScaling) const
{	double k = 1./(Gc*gScaling);
	return -k*exp(-k*delta);
}

// Get energy released (per unit volume per unit stress or Gbar/sigma) up to delta.
// gScaling must include relative toughness
double ExponentialSoftening::GetGToDelta(double delta,double gScaling) const
{	double krecip = Gc*gScaling;
	return krecip - exp(-delta/krecip)*(krecip+0.5*delta);
}

// Get G/Gc up to delta or Gbar(delta)/Gbar(infinity)
// gScaling must include relative toughness
double ExponentialSoftening::GetGoverGc(double delta,double gScaling) const
{	double k = 1./(Gc*gScaling);
	return 1. - exp(-k*delta)*(1+0.5*k*delta);
}

// dimensionless stability factor = maxSlope * sGc
double ExponentialSoftening::GetEtaStability(void) const { return 1.; }

// Find Phi(delta) function = f(delta)-delta*f'(delta))
// gscaling must include relative values
double ExponentialSoftening::GetPhiFxn(double delta,double gScaling) const
{	double k = 1./(Gc*gScaling);
	double eterm = exp(-k*delta);
	return eterm*(1+k*delta);
}

// Find Rd(delta) function = ei(f(delta)-delta*f'(delta))/(delta+ei f(delta)_^2
// gscaling and e0 must include relative values
double ExponentialSoftening::GetRdFxn(double delta,double gScaling,double e0) const
{	double k = 1./(Gc*gScaling);
	double eterm = e0*exp(-k*delta);
	double en = delta + eterm;
	return eterm*(1+k*delta)/(en*en);
}

#pragma mark ExponentialSoftening::Optional Methods

// calculate maximum cracking strain from damage (only when subcritical, initial damage state)
// gScaling and e0 must include relative values
double ExponentialSoftening::GetDeltaFromDamage(double d,double gScaling,double e0,double deltaGuess)
{	double k = 1./(Gc*gScaling);
	// to test vs. numerical code
	//#pragma omp critical (output)
	//{	cout << "# delta eq=" << gsl_sf_lambert_W0(k*d*e0/(1-d))/k << ", num="
	//	<< SofteningLaw::GetDeltaFromDamage(d,gScaling,e0,deltaGuess) << endl; }
	return gsl_sf_lambert_W0(k*d*e0/(1-d))/k;
}

// Solve for increment in crack opening strain by solving de = ddelta + e0*(f(delta+ddelta)-f(delta))
// Return ddelta if not failed or -1 if failed
// gScaling must include relative toughness, e0 must include relative strength
double ExponentialSoftening::GetDDelta(double de,double e0,double delta,double d,double gScaling) const
{   double k = 1./(Gc*gScaling);
    double x = k*e0*exp(-k*delta);
    double ddelta = de + (x + gsl_sf_lambert_W0(-x*exp(-k*de-x)))/k;
    
	// to test vs. numerical code
	//  double ddelta2 = SofteningLaw::GetDDelta(de,e0,delta,d,gScaling);
	//#pragma omp critical (output)
	//    {   cout << "# " << ddelta << "," << ddelta2 << "," << (ddelta-ddelta2) << endl;
	//    }
    
    // has it failed
    if(delta+ddelta > GetDeltaMax(gScaling)) return -1.;
    return ddelta;
}

// Find f(delta+x)/f(delta) = exp(-k x)
// gScaling must include relative toughness
double ExponentialSoftening::GetRelFFxn(double delta,double x,double gScaling) const
{	double k = 1./(Gc*gScaling);
	return exp(-k*x);
}

// Find f'(delta+x)/f(delta) = -k exp(-k x)
// gScaling must include relative toughness
double ExponentialSoftening::GetRelFpFxn(double delta,double x,double gScaling) const
{	double k = 1./(Gc*gScaling);
	return -k*exp(-k*x);
}




