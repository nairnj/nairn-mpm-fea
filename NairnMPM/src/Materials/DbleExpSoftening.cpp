/********************************************************************************
	DbleExpSoftening.hpp
	nairn-mpm-fea
	 
	Created by John Nairn, 6 July 2021.
	Copyright (c) 2021 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/DbleExpSoftening.hpp"
#include "Materials/MaterialBase.hpp"
#include "System/UnitsController.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark DbleExpSoftening::Constructors and Destructors

// Constructors
DbleExpSoftening::DbleExpSoftening() : SofteningLaw()
{
	// default to zero initial stiffness
	alpha = 0.99;				// must be > 1
	beta = 1.01;				// must be < 1
}

#pragma mark DbleExpSoftening::Initialize

// Read softening law properties
char *DbleExpSoftening::InputSofteningProperty(char *xName,int &input,double &gScaling)
{
	// Fracture toughness
	if(strcmp(xName,"alpha")==0)
	{	input=DOUBLE_NUM;
		return (char *)&alpha;
	}
	
	else if(strcmp(xName,"beta")==0)
	{	input=DOUBLE_NUM;
		return (char *)&beta;
	}
	
	return SofteningLaw::InputSofteningProperty(xName,input,gScaling);
}

// print just initiation properties to output window
void DbleExpSoftening::PrintSofteningProperties(double sigmainit)
{
	// calculate some terms
	if(alpha<=1. || beta>=1.)
	{	throw CommonException("DoubleExponential softening requires alpha>1 and beta<1",
							  "DbleExpSoftening::PrintSofteningProperties");
	}
	
	// terms used to find k (it would be 1 for single exponential)
	kterm = (alpha-beta)/(alpha*(1.-beta));
	oneOver1minusBeta = 1./(1.-beta);
	ab = alpha*beta;
	
	// ratio of peak stress to initiation stress
	if(beta>1./alpha)
	{	// peak/initiation is for delta>0
		pkRatio = (alpha-1.)/(alpha*(1.-beta)*pow(ab,1/(alpha-1.)));
	}
	else
	{	// peak is at delta=0, which is initition stress
		pkRatio = 1.;
	}
	
	// -max(f'[x]) = k*fmaxRatio
	if(beta>1/(alpha*alpha))
	{	// peak slope for delta>0
		fmaxRatio = (alpha-1.)/(alpha*(1.-beta)*pow(alpha*ab,1/(alpha-1.)));
	}
	else
	{	// max slope at delta=0
		fmaxRatio = (1.-ab)/(1.-beta);
	}
	
	// prints namem Gc, ucrit, and ecrti
	SofteningLaw::PrintSofteningProperties(sigmainit);
	
	// relative exponential decay
	MaterialBase::PrintProperty("alpha",alpha,"");
	
	// magnitude second expponential
	MaterialBase::PrintProperty("beta",beta,"");
	
	// peak stress
	double sigmapeak = sigmainit*pkRatio;
	MaterialBase::PrintProperty("s(max)",sigmapeak*UnitsController::Scaling(1.e-6),"");
	cout << endl;
}

#pragma mark DbleExpSoftening::Required Methods

// initiation law name
const char *DbleExpSoftening::GetSofteningLawName(void) const { return "Double exponential softening"; }

// This law is expnential such that f(delta) = exp(-k*delta)
// gScaling must include relative toughness
double DbleExpSoftening::GetFFxn(double delta,double gScaling) const
{	double k = kterm/(Gc*gScaling);
	return oneOver1minusBeta*(exp(-k*delta) - beta*exp(-k*alpha*delta));
}

// This law has no maximum, but effective maximum is when numerical code
// is truncated at exp(-k dmax)-beta*exp(-k*alpha*delta) = (1-beta)*minFdelta
// gScaling must include relative toughness
double DbleExpSoftening::GetDeltaMax(double gScaling) const
{	// solve exp(-x)-beta*exp(-alpha*x)-(1-beta)*min = 0
	double conTerm = (1.-beta)*minFdelta;
	double xk = -log(minFdelta);
	
	// Newton's method with g'(x) = -exp(-x)+alpha*beta*exp(-alpha*x)
	double gx = exp(-xk) - beta*exp(-alpha*xk) - conTerm;
	double gpx;
	int step = 1;
	while(step<10)
	{	// derivative derivative
		gpx = -exp(-xk) + ab*exp(-alpha*xk);
		xk -= gx/gpx;
		
		// has it converged
		gx = exp(-xk) - beta*exp(-alpha*xk) - conTerm;
		if(gx<1e-4*minFdelta) break;
		step++;
	}
	
	// x =  k*deltamax
	return xk*Gc*gScaling/kterm;
}

// This law is exponential such that such that f'(delta) = -k*exp(-k*delta)
// gScaling must include relative toughness
double DbleExpSoftening::GetFpFxn(double delta,double gScaling) const
{	double k = kterm/(Gc*gScaling);
	return -k*oneOver1minusBeta*(exp(-k*delta) - ab*exp(-k*alpha*delta));
}

// Get energy released (per unit volume per unit stress or Gbar/sigma) up to delta.
// gScaling must include relative toughness
double DbleExpSoftening::GetGToDelta(double delta,double gScaling) const
{	double krecip = Gc*gScaling;
	double k = kterm/krecip;
	double ka = k*alpha;
	double hd = 0.5*delta;
	return krecip - oneOver1minusBeta*(exp(-k*delta)*(1./k+hd) - beta*exp(-ka*delta)*(1./ka+hd));
}

// Get G/Gc up to delta or Gbar(delta)/Gbar(infinity)
// gScaling must include relative toughness
double DbleExpSoftening::GetGoverGc(double delta,double gScaling) const
{	double k = kterm/(Gc*gScaling);
	double ka = k*alpha;
	double hd = 0.5*delta;
	return 1. - (ka/(alpha-beta))*(exp(-k*delta)*(1./k+hd) - beta*exp(-ka*delta)*(1./ka+hd));
}

// dimensionless stability factor = 1/(maxSlope * sGc)
double DbleExpSoftening::GetEtaStability(void) const
{	return 1./(kterm*fmaxRatio);
}

// Find Phi(delta) function = f(delta)-delta*f'(delta))
// gscaling must include relative values
double DbleExpSoftening::GetPhiFxn(double delta,double gScaling) const
{	double k = kterm/(Gc*gScaling);
	double ka = k*alpha;
	return oneOver1minusBeta*(exp(-k*delta)*(1.+k*delta) - beta*exp(-ka*delta)*(1.+ka*delta));
}

// Find Rd(delta) function = ei(f(delta)-delta*f'(delta))/(delta+ei f(delta)_^2
// gscaling and e0 must include relative values
double DbleExpSoftening::GetRdFxn(double delta,double gScaling,double e0) const
{	double k = kterm/(Gc*gScaling);
	double ka = k*alpha;
	double ed = e0*exp(-k*delta);
	double ekd = e0*exp(-ka*delta);
	double en = delta + oneOver1minusBeta*(ed - beta*ekd);
	return oneOver1minusBeta*(ed*(1.+k*delta)-beta*ekd*(1.+ka*delta))/(en*en);
}
