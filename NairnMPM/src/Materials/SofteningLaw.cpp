/********************************************************************************
	SofteningLaw.cpp
	nairn-mpm-fea

	Abstract base class for softening laws
 
	Created by John Nairn, June 26, 2015.
	Copyright (c) 2008 John A. Nairn, All rights reserved.
*******************************************************************************/

#include "stdafx.h"
#include "Materials/SofteningLaw.hpp"
#include "Materials/MaterialBase.hpp"
#include "System/UnitsController.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Exceptions/CommonException.hpp"

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
    unscaledDeltaMax = -log(minFdelta)*Gc;			// only used by exponential softening
	double ucrit = GetDeltaMax(1.)/sigmac;
	MaterialBase::PrintProperty("ucrit",ucrit,UnitsController::Label(CULENGTH_UNITS));
	
	// critical strain < ucrit/(min particle length)
	double ecrit = ucrit/mpmgrid.GetGlobalMinParticleLength();
	MaterialBase::PrintProperty("ecrit>",100.*ecrit,"%");
	cout << endl;
}

#pragma mark SofteningLaw::General Methods (optional overrides)

// Get energy released (per unit volume per unit stress or Gbar/sigma) up to delta.
// General result int_0^delta f(delta) - delta*f(delta)/2
// Current code not longer calls this method. Cohesive laws need not override
double SofteningLaw::GetGToDelta(double delta,double gScaling) const
{	throw CommonException("GetGToDelta() called for cohesive law lacking its implementation",
						  "SofteningLaw::GetGoverGc");
}

// Get G/Gc up to delta Gbar(delta)/Gbar(deltaMax) = 0.5*delta/(0.5*deltaMax) = delta/deltaMax
// Also equal to GetGToDelta(delta)/GetGToDelta(deltaMax)
// gScaling must include relative toughness
// This method is only called for cubic failure surface and for isotropic materials, it is only called in 3D
// Laws to be used in these modes must override.
double SofteningLaw::GetGoverGc(double delta,double gScaling) const
{	throw CommonException("GetGOverGc() called for cohesive law lacking its implementation",
						  "SofteningLaw::GetGoverGc");
}

// Find f(delta+x)/f(delta)
// Only needed in Newton's law in GetDDelta() and law can override if has more efficient answer
// 		but need not override is not making use of Newton's method
// gScaling must include relative toughness, e0 must include relative strength
double SofteningLaw::GetRelFFxn(double delta,double x,double gScaling) const
{	return GetFFxn(delta+x,gScaling)/GetFFxn(delta,gScaling);
}

// Find f'(delta+x)/f(delta)
// Only needed in Newton's law in GetDDelta() and law can override if has more efficient answer
// 		but need not override is not making use of Newton's method
// gScaling must include relative toughness
double SofteningLaw::GetRelFpFxn(double delta,double x,double gScaling) const
{	return GetFpFxn(delta+x,gScaling)/GetFFxn(delta,gScaling);
}

// Calculate delta parameter from damage parameter D  (only when subcritical, initial damage state)
// numerically solves g(delta) = 0 =  delta*(1-d) - d*e0*f(delta)
// For Newton's method g'(delta) = 1-d - d*e0*f'(delta)
// g(delta+ddelta) = 0 = g(delta) + g'(delta)*ddelta or ddelta = -g(delta)/g'(delta)
// gScaling and e0 must include relative values
double SofteningLaw::GetDeltaFromDamage(double d,double gScaling,double e0,double deltaGuess)
{	// if outside brackets, return and end point
	double dmin = 0.;
	if(d<=0.) return dmin;
	double dmax = GetDeltaMax(gScaling);
	if(d>=1.) return dmax;
	
	// solution for delta bracketed between 0 (g(0) = -d*e0 < 0) and
	//			deltaMax (g(deltaMax) = delta*(1-d) > 0)
	
	/*
	// Simple binary search, Newton's method might be better, but this only
	// called during initiatialization phaase
	double dmid;
	int i=0,nsteps = 10;
	while(i<nsteps)
	{	// value at midpoint
		dmid = 0.5*(dmin+dmax);
		double gmid = d*e0*GetFFxn(dmid,gScaling) - dmid*(1.-d);
		
		// see if done or move one bracket
		if(gmid==0.)
			return dmid;
		else if(gmid>0.)
			dmin = dmid;
		else
			dmax = dmid;
		
		// next step
		i++;
	}
	
	// last step - interpolate g = gmin + (delta-dmin)*(gmax-gmin)/(dmax-dmin) = gmin + fract*(delta-dmin)
	// Sovling for zero: delta = dmin - gmin/fract = (dmin*gmax-gmin*dmax)/(gmax-gmin)
	double gmin = d*e0*GetFFxn(dmin,gScaling) - dmin*(1.-d);
	double gmax = d*e0*GetFFxn(dmax,gScaling) - dmax*(1.-d);
	dmid = (dmin*gmax-gmin*dmax)/(gmax-gmin);
	return dmid;
	*/
	
	// This is unbracketed Newton with no checks for non-convergence. It is likely
	// safe for typical softening laws. It is not used for linear (or any other model
	// with an analytical solution for delta as function of d)
	double delk = deltaGuess>-0. ? deltaGuess : 0.5*dmax;
	double dx;
	int step = 1;
	while(true)
	{	// update iterative variables (lambda, alpha)
		double gdel = delk*(1.-d) - d*e0*GetFFxn(delk,gScaling);
		double slope = 1. - d - d*e0*GetFpFxn(delk,gScaling);
		
		dx = gdel/slope;
		double temp = delk;
		delk -= dx;
		if(temp == delk) break;  // change in root is negligible
		
		// check for convergence
		step++;
		if(step>15 || fabs(dx/dmax)<1.e-8) break;
	}

	return delk;
}

//#define TRACK_NEWTON

// Solve for increment in crack opening strain by solving de = ddelta + e0*(f(delta+ddelta)-f(delta))
//      for ddelta
// Solve numerically here by Newton's method with custom bracketing.
//      Subclass can override if has better option
// Return ddelta if not failed or -1 if failed
double SofteningLaw::GetDDelta(double de,double e0,double delta,double d,double gScaling) const
{
    // The solution for ddelta is bracketd by D*de (because f'(delta)<f(delta)/delta) <= 0 )
    // and de + e0*f(delta) (if ddelta caused failure where f(delta+ddelta) = 0).
	double xl = d*de;
	double fdelta = GetFFxn(delta,gScaling);
	double b = e0*fdelta;
	double xh = de + b;
	double bdenom = 1./b;
	
	// initial guess in the middle
	int step = 1;
	double xk = xl + 0.05*b;
	double gpx,dx,range,fx = GetRelFFxn(delta,xk,gScaling);
	double gx = bdenom*(xk - de) + fx - 1.;
	
	while(true)
	{	// update iterative variables (lambda, alpha)
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
		fx = GetRelFFxn(delta,xk,gScaling);		// final f(delta+xk)/f(delta)
		gx = bdenom*(xk - de) + fx - 1.;
#ifdef TRACK_NEWTON
		cout << gx ;
#endif
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

//#define TRACK_NEWTON

// Solve for increment in delta during elastic deformation based on starting d
//		The equation is delta + ddeltae = F(delta + ddeltae,alpa+dalpha)/(E(1-d))
// On input scaleAlpha is updated value and sigmaAlpha is updated stress.
// Only needed when softening law depends on extra variables
// Solve numerically by Newton's method with bracketing; subclass can override if better solution
double SofteningLaw::GetDDeltaElastic(double delta,double sigmaAlpha,double scaleAlpha,double d,double E1md) const
{
    // negative stress means the softening surface does not depend on other variables
    // delta=0 has trivial solution of zero
    if(sigmaAlpha<0. || delta==0.) return 0.;
    
    // solving g(x) = E1md*(delta+x) - d*sigmaAlpha*GetFFxn(delta+x,scaleAlpha) = 0
    //    with
    // g'(x) = E1md - d*sigmaAlpha*GetFpFxn(delta+x,scaleAlpha)
    
    // Solution bracketed by xl <= x <= xh with g(xl)<0 and g(xh)>0
    double xl,xh;
    double diff = E1md*delta - d*sigmaAlpha*GetFFxn(delta,scaleAlpha);
    double deltaMax = GetDeltaMax(scaleAlpha);
    if(diff==0.)
    {   // for dAlpha=0 or for law independent of alpha
        return 0.;
    }
    else if(diff<0.)
    {   xl = 0.;
        xh = deltaMax - delta;
    }
    else
    {   xl = -delta;
        xh = 0.;
    }
    double deltaConverge = 1.e-9*deltaMax;
    
    // initial guess in the middle
    int step = 1;
    double xk = 0.,gpx,dx,dxold=fabs(xh-xl);
    double gx = E1md*(xk+delta) - d*sigmaAlpha*GetFFxn(delta+xk,scaleAlpha);

#ifdef TRACK_NEWTON
    cout << "#... " << xl << " to " << xh << " with " << gx << " at " << xk << endl;
#endif

    while(true)
    {   // update derivative
        gpx = E1md - d*sigmaAlpha*GetFpFxn(delta+xk,scaleAlpha);
        
#ifdef TRACK_NEWTON
        cout << "#   " << step << ": gx = " << gx << " : ";
#endif
        //if( (((xk-xl)*gpx-gx)*((xk-xh)*gpx-gx) >= 0.) || (fabs(2.*gx)>fabs(dxold*gpx)) )
        if( (((xk-xl)*gpx-gx)*((xk-xh)*gpx-gx) >= 0.) )
        {   // means this jump is not between xl and xh
            
            // seen to occur when get close enough to solution, so accept it
            if(DbleEqual(gx,0.)) break;
            
            // take midpoint instead
            dxold = dx;
            dx = 0.5*(xh-xl);
            xk = xl+dx;
#ifdef TRACK_NEWTON
            cout << "B (" << xl << "," << xk << "," << xh << ") (" <<
                ((xk-xl)*gpx-gx)*((xk-xh)*gpx-gx) << " or " << (fabs(2.*gx)>fabs(dxold*gpx)) << ") ";
#endif
            if(xl == xk) break;    // change in root is negligible
        }
        else
        {   dxold = dx;
            double xkprev = xk;
            dx = gx/gpx;
            xk -= dx;
#ifdef TRACK_NEWTON
            cout << "N (" << xl << "," << xk << "," << xh << ") ";
#endif
            if(xk == xkprev) break;  // change in root is negligible
        }
        
        // is it converged
        if(fabs(dx)<deltaConverge || step>15) break;
        
        // next step: new limits and increment step
        gx = E1md*(xk+delta) - d*sigmaAlpha*GetFFxn(delta+xk,scaleAlpha);
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
    cout << " done: xk=" << xk << ", gx=" << (E1md*(xk+delta) - d*sigmaAlpha*GetFFxn(delta+xk,scaleAlpha)) << endl;
#endif
    
    return xk;
}

#pragma mark SofteningLaw::General Accessors (Optional overrides)

// toughness
double SofteningLaw::GetGc(void) const { return Gc; }

// if it linear softening (constant derivative)
bool SofteningLaw::IsLinear(void) const { return false; }

// Check if new delta is too high (adjust ddel if it is)
// gScaling must include relative toughness
bool SofteningLaw::HasFailed(double delta,double &ddel,double gScaling) const
{	double deltaMax = GetDeltaMax(gScaling);
	if(delta+ddel>deltaMax)
	{	ddel=deltaMax-delta;
		return true;
	}
	return false;
}

