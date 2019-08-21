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

#pragma mark ExponentialSoftening::Methods

// This law is expnential such that f(delta) = exp(-k*delta)
double ExponentialSoftening::GetFFxn(double delta,double gScaling) const
{	double k = 1./(Gc*gScaling);
	return exp(-k*delta);
}

// This law is exponential such that such that f'(delta) = -k*exp(-k*delta)
double ExponentialSoftening::GetFpFxn(double delta,double gScaling) const
{	double k = 1./(Gc*gScaling);
	return -k*exp(-k*delta);
}

// Find f(delta+x)/f(delta) = exp(-k x)
double ExponentialSoftening::GetRelFFxn(double delta,double x,double gScaling) const
{	double k = 1./(Gc*gScaling);
	return exp(-k*x);
}

// Find f'(delta+x)/f(delta) = -k exp(-k x)
double ExponentialSoftening::GetRelFpFxn(double delta,double x,double gScaling) const
{	double k = 1./(Gc*gScaling);
	return -k*exp(-k*x);
}

// Get energy released (per unit volume per unit stress or Gbar/sigma) up to delta.
double ExponentialSoftening::GetGToDelta(double delta,double gScaling) const
{	double krecip = Gc*gScaling;
	return krecip - exp(-delta/krecip)*(krecip+0.5*delta);
}

// Get G/Gc up to delta or Gbar(delta)/Gbar(infinity)
double ExponentialSoftening::GetGoverGc(double delta,double gScaling) const
{	double k = 1./(Gc*gScaling);
	return 1. - exp(-k*delta)*(1+0.5*k*delta);
}

// Get maximium decreasing slope (max(-f'(delta)) for this softening law
double ExponentialSoftening::GetMaxSlope(double gScaling) const
{	return 1./(Gc*gScaling);
}

#pragma mark ExponentialSoftening::Accessors

// This law has no maximum, but effective maximum is when numerical code
// is truncated at exp(-k dmax) = minFdelta or dmax = - ln(minFdelta)/k
double ExponentialSoftening::GetDeltaMax(double gScaling) const
{	return -log(minFdelta)*Gc*gScaling;
}

// initiation law name
const char *ExponentialSoftening::GetSofteningLawName(void) const { return "Exponential softening"; }
