/********************************************************************************
	LinearSoftening.cpp
	nairn-mpm-fea

	Softening law is linear f(delta) = 1 - delta/deltaMax

	but deltaMax depends on particle size and crack area by

 		deltaMax = (2Ac/(Vp rho sigma(sp)) * Gc
 		umax = Vp deltaMax/Ac = 2 Gc/(rho sigma(sp))

	This scaling factor is input to functions that need deltaMax

 		gscaling = Ac/(Vp rho sigma(sp))
 		deltaMax = 2*Gc*gscaling

	Created by John Nairn, June 26, 2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.
*******************************************************************************/

#include "stdafx.h"
#include "Materials/LinearSoftening.hpp"
#include "Materials/MaterialBase.hpp"
#include "System/UnitsController.hpp"

#pragma mark LinearSoftening::Methods

// Linear law is f(delta) = 1 - delta/deltaMax
// gScaling must include relative toughness
double LinearSoftening::GetFFxn(double delta,double gScaling) const
{	double deltaMax = 2.*Gc*gScaling;
	return fmax(1.-delta/deltaMax,0.);
}

// Linear derivative is constant f'(delta) = -1/deltaMax
// gScaling must include relative toughness
double LinearSoftening::GetFpFxn(double delta,double gScaling) const
{	double deltaMax = 2.*Gc*gScaling;
	return delta<=deltaMax ? -1./deltaMax : 0. ;
}

// Solve for increment in crack opening strain by solving de = ddelta + e0*(f(delta+ddelta)-f(delta))
// Return ddelta if not failed or -1 if failed
// gScaling must include relative toughness, e0 must include relative strength
double LinearSoftening::GetDDelta(double de,double e0,double delta,double gScaling) const
{	double deltaMax = 2.*Gc*gScaling;
	double ddelta = de/(1.-e0/deltaMax);
	
	// has it failed
	if(delta+ddelta>deltaMax) ddelta = -1.;
	return ddelta;
}

// Solve for increment in delta during elastic deformation based on starting d
// Only needed when softening law depends on extra variables
// Solve numerically by Newton's method with bracketing; subclass can override if better solution
double LinearSoftening::GetDDeltaElastic(double delta,double sigmaAlpha,double scaleAlpha,double d,double E1md) const
{
    // negative stress means the softening surface does not depend on other variables
    // delta=0 has trivial solution of zero
    if(sigmaAlpha<0. || delta==0.) return 0.;
    double deltaMax = 2.*Gc*scaleAlpha;
    double dsigma = d*sigmaAlpha;
    //cout << "# Eq=" << (dsigma*deltaMax/(E1md*deltaMax+dsigma) - delta) << " vs "
    //            << SofteningLaw::GetDDeltaElastic(delta,sigmaAlpha,scaleAlpha,d,E1md) << endl;
    return dsigma*deltaMax/(E1md*deltaMax+dsigma) - delta;
}

// Get energy released (per unit volume per unit stress or Gbar/sigma) up to delta.
// General result int_0^delta f(delta) - delta*f(delta)/2
// For linear softening (delta - delta^2/(2*deltaMax) - delta*fval/2 where fval = 1-delta/deltaMax
//     = delta*(1 - delta/(2*deltaMax)) - delta*(1/2 - delta/(2*deltaMax)) = delta/2
double LinearSoftening::GetGToDelta(double delta,double gScaling) const
{	return 0.5*delta;
}

// Get G/Gc up to delta Gbar(delta)/Gbar(deltaMax) = 0.5*delta/(0.5*deltaMax) = delta/deltaMax
// Also equal to GetGToDelta(delta)/GetGToDelta(deltaMax)
// gScaling must include relative toughness
double LinearSoftening::GetGoverGc(double delta,double gScaling) const
{	double deltaMax = 2.*Gc*gScaling;
	return delta/deltaMax;
}

// Get maximium decreasing slope (max(-f'(delta)) for this softening law
// Linear is 1/deltamax
// gScaling must include relative toughness
double LinearSoftening::GetMaxSlope(double gScaling) const
{	double deltaMax = 2.*Gc*gScaling;
	return 1./deltaMax;
}

// calculate maximum cracking strain from damage (only when subcritical, initial damage state)
// gScaling and e0 must include relative values
double LinearSoftening::GetDeltaFromDamage(double d,double gScaling,double e0,double deltaGuess)
{	double deltaMax = 2.*Gc*gScaling;
	
	// to test numerical code
	//cout << "# delta eq=" << (deltaMax*d*e0/(deltaMax*(1-d) + e0*d)) << ", num="
	//		<< SofteningLaw::GetDeltaFromDamage(d,gScaling,e0,deltaGuess) << endl;
	
	return deltaMax*d*e0/(deltaMax*(1-d) + e0*d);
}

// Find Rd(delta) function = ei(f(delta)-delta*f'(delta))/(delta+ei f(delta)_^2
// gscaling and e0 must include relative values
double LinearSoftening::GetRdFxn(double delta,double gScaling,double e0) const
{	double deltaMax = 2.*Gc*gScaling;
	double en = delta + e0*(1.-delta/deltaMax);
	return e0/(en*en);
}

// Find Phi(delta) function = f(delta)-delta*f'(delta))
// gscaling must include relative values
double LinearSoftening::GetPhiFxn(double delta,double gScaling) const
{	return 1.;
}

#pragma mark LinearSoftening::Accessors

// maximum delta
// gScaling must include relative toughness
double LinearSoftening::GetDeltaMax(double gScaling) const { return 2.*Gc*gScaling; }

// initiation law name
const char *LinearSoftening::GetSofteningLawName(void) const { return "Linear softening"; }

// if it linear softening (constant derivative)
bool LinearSoftening::IsLinear(void) const { return true; }

// dimensionless stability factor
double LinearSoftening::GetEtaStability(void) const { return 2.; }


