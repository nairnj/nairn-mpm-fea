/********************************************************************************
	LinearSoftening.cpp
	nairn-mpm-fea

	Softening law is linear f(delta) = 1 - delta/deltaMax

	but deltaMax depends on particle size and crack area by

	deltaMax = (2Ac/(Vp rho sigma(sp)) * Gc
	umax = Vp deltaMax/Ac = 2 Gc/(rho sigma(sp))

	This scaling factor is input to function that need deltaMax

	gscaling = Ac/(Vp rho sigma(sp))
	deltaMax = 2*Gc*gscaling

	Created by John Nairn, June 26, 2015 (transferred to this code Jan 21, 2017)
	Copyright (c) 2015 John A. Nairn, All rights reserved.
*******************************************************************************/

#include "stdafx.h"
#include "Materials/LinearSoftening.hpp"
#include "Materials/MaterialBase.hpp"
#include "System/UnitsController.hpp"

#pragma mark LinearSoftening::Methods

// Linear law is f(delta) = 1 - delta/deltaMax
double LinearSoftening::GetFFxn(double delta,double gScaling) const
{	double deltaMax = 2.*Gc*gScaling;
	return fmax(1.-delta/deltaMax,0.);
}

// Linear derivative is constant f'(delta) = -1/deltaMax
double LinearSoftening::GetFpFxn(double delta,double gScaling) const
{	double deltaMax = 2.*Gc*gScaling;
	return delta<=deltaMax ? -1./deltaMax : 0. ;
}

// Solve for increment in crack opening strain by solving de = ddelta + e0*(f(delta+ddelta)-f(delta))
// Return ddelta if not failed or -1 if failed
double LinearSoftening::GetDDelta(double de,double e0,double delta,double gScaling) const
{	double deltaMax = 2.*Gc*gScaling;
	double ddelta = de/(1.-e0/deltaMax);
	
	// has it failed
	if(delta+ddelta>deltaMax) ddelta = -1.;
	return ddelta;
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
double LinearSoftening::GetGoverGc(double delta,double gScaling) const
{	double deltaMax = 2.*Gc*gScaling;
	return delta/deltaMax;
}

// Get maximium decreasing slope (max(-f'(delta)) for this softening law
// Linear is 1/deltamax
double LinearSoftening::GetMaxSlope(double gScaling) const
{	double deltaMax = 2.*Gc*gScaling;
	return 1./deltaMax;
}

// calculate maximum cracking strain from damage
double LinearSoftening::GetDeltaFromDamage(double d,double gScaling,double e0)
{	double deltaMax = 2.*Gc*gScaling;
	
	// to test numerical code
	//cout << "# " << (deltaMax*d*e0/(deltaMax*(1-d) + e0*d)) << ", " << SofteningLaw::GetDeltaFromDamage(d,gScaling,e0) << endl;
	
	return deltaMax*d*e0/(deltaMax*(1-d) + e0*d);
}

#pragma mark LinearSoftening::Accessors

// maximum delta
double LinearSoftening::GetDeltaMax(double gScaling) const { return 2.*Gc*gScaling; }

// initiation law name
const char *LinearSoftening::GetSofteningLawName(void) const { return "Linear softening"; }

// if it linear softening (constant derivative)
bool LinearSoftening::IsLinear(void) const { return true; }

