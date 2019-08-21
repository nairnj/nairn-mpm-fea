/********************************************************************************
	SofteningLaw.cpp
	nairn-mpm-fea

	Softening law is linear f(delta) = 1 - delta/deltaMax
 
	but deltaMax depends on particle size and crack area by
 
		deltaMax = (2Ac/(Vp rho sigma(sp)) * Gc
		umax = Vp deltaMax/Ac = 2 Gc/(rho sigma(sp))
 
	This scaling factor is input to function that need deltaMax

		gscaling = Ac/(Vp rho sigma(sp))
		deltaMax = 2*Gc*gscaling
 
	Created by John Nairn, June 26, 2015.
	Copyright (c) 2008 John A. Nairn, All rights reserved.
*******************************************************************************/

#include "stdafx.h"
#include "Materials/SofteningLaw.hpp"
#include "Materials/MaterialBase.hpp"
#include "System/UnitsController.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"

// Tested with some laws and wide range seems to work
// Likely no need to make it smaller
#define NORM_GX_CONVERGE 1.e-6

#pragma mark SofteningLaw::Constructors and Destructors

// Constructors
SofteningLaw::SofteningLaw()
{
	Gc = 1.e40;

	// To avoid numerical issues in numerical options
	minFdelta = 0.01;
}

// Destructor (and it is virtual)
SofteningLaw::~SofteningLaw() {}


#pragma mark SofteningLaw::Initialize

// Read softening law properties
char *SofteningLaw::InputSofteningProperty(char *xName,int &input,double &gScaling)
{
	// Fracture toughness
    if(strcmp(xName,"Gc")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&Gc,gScaling,1000.);
	}

	// minimum F
	else if(strcmp(xName,"min")==0)
	{	input=DOUBLE_NUM;
		return (char *)&minFdelta;
	}

    // is not a softening law property
    return NULL;
}

// print just initiation properties to output window
void SofteningLaw::PrintSofteningProperties(double sigmac)
{
	cout << GetSofteningLawName() << endl;
	MaterialBase::PrintProperty("Gc",Gc*UnitsController::Scaling(1.e-3),UnitsController::Label(ERR_UNITS));

	// ucrit = d(max)/(s sigmac) and we assume d(max) is proportional to s
	double ucrit = GetDeltaMax(1.)/sigmac;
	MaterialBase::PrintProperty("ucrit",ucrit,UnitsController::Label(CULENGTH_UNITS));
	
	// critical strain < ucrit/(min particle length)
	double ecrit = ucrit/mpmgrid.GetGlobalMinParticleLength();
	MaterialBase::PrintProperty("ecrit>",100.*ecrit,"%");
	cout << endl;
}

#pragma mark SofteningLaw::Methods

// Find f(delta+x)/f(delta)
// Only needed in Newton's law in GetDDelta() and law can override if has more efficient answer
double SofteningLaw::GetRelFFxn(double delta,double x,double gScaling) const
{	return GetFFxn(delta+x,gScaling)/GetFFxn(delta,gScaling);
}

// Find f'(delta+x)/f(delta)
// Only needed in Newton's law in GetDDelta() and law can override if has more efficient answer
double SofteningLaw::GetRelFpFxn(double delta,double x,double gScaling) const
{	return GetFpFxn(delta+x,gScaling)/GetFFxn(delta,gScaling);
}

// Check if new delta is too high (adjust ddel if it is)
bool SofteningLaw::HasFailed(double delta,double &ddel,double gScaling) const
{	double deltaMax = GetDeltaMax(gScaling);
	if(delta+ddel>deltaMax)
	{	ddel=deltaMax-delta;
		return true;
	}
	return false;
}

//#define TRACK_NEWTON

// Solve for increment in crack opening strain by solving de = ddelta + e0*(f(delta+ddelta)-f(delta))
//      for ddelta
// Solve numerically here by Newton's method with custom bracketing.
//      Subclass can override if has better option
// Return ddelta if not failed or -1 if failed
double SofteningLaw::GetDDelta(double de,double e0,double delta,double gScaling) const
{
    // The solution for ddelta is bracketd by de (because e0*(f(delta+ddelta)-f(delta)) <= 0 )
    // and de + e0*f(delta) (if ddelta caused failure where f(delta+ddelta) = 0).
	double xl = de;
	double fdelta = GetFFxn(delta,gScaling);
	double b = e0*fdelta;
	double xh = de + b;
	double bdenom = 1./b;
	
	// initial guess in the middle
	int step = 1;
	double xk = xl + 0.05*b;
	double gx, gpx, dx, range, fx = GetRelFFxn(delta,xk,gScaling);
	
	while(true)
	{	// update iterative variables (lambda, alpha)
		gx = bdenom*(xk - de) + fx - 1.;
		gpx = bdenom + GetRelFpFxn(delta,xk,gScaling);
		
		// We assume softening law is well behaved and therefore
		//    extending outside range is round off error. When happens
		//    set to new value just inside the range
		// Proposed new value is xnew = xk - gx/gpx
		// If below xl then xnew-xl<0 or (xk-xl)*gpx-gx<0
		// If above xh then xnew-xh>0 or (xk-xh)*gpx-gx>0
#ifdef TRACK_NEWTON
		if(step==1)
			cout << "#..." << step << ": ";
		else
			cout << "#   " << step << ": ";
#endif
		range = xh - xl;
		if(((xk-xl)*gpx-gx) < 0.)
		{	xk = xl + 1.e-6*range;
#ifdef TRACK_NEWTON
			cout << "BL (" << xl << "," << xk << "," << xh << ") ";
#endif
			if(xl == xk) break;    // change in root is negligible
		}
		else if(((xk-xh)*gpx-gx) > 0.)
		{	xk = xh - 1.e-6*range;
#ifdef TRACK_NEWTON
			cout << "BH (" << xl << "," << xk << "," << xh << ") ";
#endif
			if(xh == xk) break;    // change in root is negligible
		}
		else
		{	double xkprev = xk;
			dx = gx/gpx;
			xk -= dx;
#ifdef TRACK_NEWTON
			cout << "N ";
#endif
			if(xk == xkprev) break;  // change in root is negligible
		}
		
		// is it converged
#ifdef TRACK_NEWTON
		cout << gx ;
#endif
		fx = GetRelFFxn(delta,xk,gScaling);		// final f(delta+xk)/f(delta)
		if(fabs(gx)<NORM_GX_CONVERGE || step>15) break;
		
		// next step: new limits and increment step
		if(gx < 0.)
			xl = xk;
		else
			xh = xk;
		step++;
#ifdef TRACK_NEWTON
		cout << endl;
#endif
	}

#ifdef TRACK_NEWTON
	cout << " done: " << xk << "," << gx << "," << fx*fdelta << endl;
#endif
	
	// failed if fx is close enough to zero
	if(fx*fdelta<minFdelta) return -1.;
	return xk;
}

// calculate maximum cracking strain from damage
// numerically solves g(delta) = 0 = d*e0*f(delta) - delta*(1-d)
// For Newton's method g'(delta) = d*e0*f'(delta) - 1 + d or increment ddelta = -g(delta)/(d*e0*f'(delta) - 1 + d)
// g(delta+ddelta) = 0 = g(delta) + g'(delta)*ddelta or ddelta = -g(delta)/g'(delta)
double SofteningLaw::GetDeltaFromDamage(double d,double gScaling,double e0)
{	double dmin = 0.;
	if(d<=0.) return dmin;
	double dmax = GetDeltaMax(gScaling);
	if(d>=1.) return dmax;
	
	// solution for delta bracketed between 0 (g(0) = d*e0 >0) and
	//			deltaMax (g(deltaMax) = -delta*(1-d) < 0)
	
	// Simple binary search, Newton's method might be better, but this only
	// called during initiatialization phaase
	
	int i=0,nsteps = 10;
	while(i<nsteps)
	{	// value at midpoint
		double dmid = 0.5*(dmin+dmax);
		double gmid = d*e0*GetFFxn(dmid,gScaling) - dmid*(1.-d);
		
		// see if done or move one bracket
		if(gmid==0.)
			return dmid;
		else if(gmid>0.)
			dmin = dmid;
		else
			dmax = dmid;
		
		// next step
		nsteps++;
	}
	
	// last step - interpolate g = gmin + (delta-dmin)*(gmax-gmin)/(dmax-dmin) = gmin + fract*(delta-dmin)
	// Sovling for zero: delta = dmin - gmin/fract = (dmin*gmax-gmin*dmax)/(gmax-gmin)
	double gmin = d*e0*GetFFxn(dmin,gScaling) - dmin*(1.-d);
	double gmax = d*e0*GetFFxn(dmax,gScaling) - dmax*(1.-d);
	return (dmin*gmax-gmin*dmax)/(gmax-gmin);
}

#pragma mark SofteningLaw::Accessors

// toughness
double SofteningLaw::GetGc(void) const { return Gc; }

// if it linear softening (constant derivative)
bool SofteningLaw::IsLinear(void) const { return false; }

