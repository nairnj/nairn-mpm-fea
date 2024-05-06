/********************************************************************************
	Nonlinear2Hardening.hpp
	nairn-mpm-fea

	Created by John Nairn, 2/8/2103
	Copyright (c) 2013 John A. Nairn, All rights reserved.

	Hardening Law is
		sigma = yield (1 + beta ep^npow)
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Materials/Nonlinear2Hardening.hpp"

#pragma mark Nonlinear2Hardening::Constructors and Destructors

// Constructors
Nonlinear2Hardening::Nonlinear2Hardening() {}

// Constructors
Nonlinear2Hardening::Nonlinear2Hardening(MaterialBase *pair) : NonlinearHardening(pair)
{
	lawID = NONLINEAR2HARENDING_ID;
}

// get reduced stress than done
const char *Nonlinear2Hardening::VerifyAndLoadProperties(int np)
{
	// call first to get reduced yield stress
	HardeningLawBase::VerifyAndLoadProperties(np);
	
	// maximum alpha when softening
	if(beta<0.)
		alphaMax = pow((yldredMin/yldred - 1.)/beta , 1./npow);
	
	// base call above never has an error
	return NULL;
}
#pragma mark NonlinearHardening::Law Methods

// Return yield stress for current conditions (alpint for cum. plastic strain and dalpha/delTime for plastic strain rate)
double Nonlinear2Hardening::GetYield(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	return a->alpint < alphaMax ? yldred*(1.+beta*pow(a->alpint,npow)) : yldredMin ;
}

// Get derivative of sqrt(2./3.)*yield with respect to lambda for plane strain and 3D
// ... and using dep/dlambda = sqrt(2./3.)
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)lamda or depdot/dlambda = sqrt(2./3.)/delTime
double Nonlinear2Hardening::GetKPrime(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	return a->alpint < alphaMax ? TWOTHIRDS*yldred*beta*npow*pow(a->alpint,npow-1.) : 0. ;
}

// Get derivative of (1./3.)*yield^2 with respect to lambda for plane stress only
// ... and using dep/dlambda = sqrt(2./3.)*fnp1 where ep=alpint
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)*lambda*fnp1 or depdot/dlambda = sqrt(2./3.)*fnp1/delTime
// Also equal to sqrt(2./3.)*GetYield()*GetKPrime()*fnp1, but in separate call for efficiency
double Nonlinear2Hardening::GetK2Prime(MPMBase *mptr,double fnp1,double delTime,HardeningAlpha *a,void *properties) const
{
    if(DbleEqual(a->alpint,0.)) return 0.;
	if(a->alpint < alphaMax)
	{	double alphan = pow(a->alpint,npow);
		return SQRT_EIGHT27THS*yldred*yldred*beta*npow*(1.+beta*alphan)*alphan*fnp1/a->alpint;
	}
	else
		return 0.;
}

#pragma mark NonlinearHardening::Accessors

// hardening law name
const char *Nonlinear2Hardening::GetHardeningLawName(void) const { return "Nonlinear hardening 2"; }

