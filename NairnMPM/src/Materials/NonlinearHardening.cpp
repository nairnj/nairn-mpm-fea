/********************************************************************************
    NonlinearHardening.hpp
    nairn-mpm-fea

    Created by John Nairn, 1/17/2103
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Hardening Law is
        sigma = yield (1 + beta ep)^npow
********************************************************************************/

#include "stdafx.h"
#include "Materials/NonlinearHardening.hpp"
#include "System/UnitsController.hpp"

#pragma mark NonlinearHardening::Constructors and Destructors

// Constructors
NonlinearHardening::NonlinearHardening() {}

// Constructors
NonlinearHardening::NonlinearHardening(MaterialBase *pair) : HardeningLawBase(pair)
{
	beta = 0.;
	npow = 1.;
	alphaMax = 1.e50;
}

#pragma mark LinearHardening::Initialize

// Read hardening law properties
char *NonlinearHardening::InputHardeningProperty(char *xName,int &input,double &gScaling)
{
    // Khard: coefficient of plastic strains for non-linear hardening (beta)
    if(strcmp(xName,"Khard")==0)
    {   input=DOUBLE_NUM;
        return((char *)&beta);
    }
	
    // mhard: power in non-linear hardening
    else if(strcmp(xName,"nhard")==0)
    {   input=DOUBLE_NUM;
        return((char *)&npow);
    }
    
    return HardeningLawBase::InputHardeningProperty(xName,input,gScaling);
}

// get reduced stress than done
const char *NonlinearHardening::VerifyAndLoadProperties(int np)
{
	// call first to get reduced yield stress
	HardeningLawBase::VerifyAndLoadProperties(np);
	
	// maximum alpha when softening
	if(beta<0.)
		alphaMax = (pow(yldredMin/yldred, 1./npow) - 1.)/beta;
	
	// base call above never has an error
	return NULL;
}

// print just yield properties to output window
void NonlinearHardening::PrintYieldProperties(void) const
{
    cout << GetHardeningLawName() << endl;
    MaterialBase::PrintProperty("yld",yield*UnitsController::Scaling(1.e-6),"");
    MaterialBase::PrintProperty("beta",beta,"");
    MaterialBase::PrintProperty("n",npow,"");
	if(beta<0.)
		MaterialBase::PrintProperty("yldMin",yieldMin*UnitsController::Scaling(1.e-6),"");
    cout << endl;
}

#pragma mark NonlinearHardening::Law Methods

// Return yield stress for current conditions (alpint for cum. plastic strain and dalpha/delTime for plastic strain rate)
double NonlinearHardening::GetYield(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{   
	return a->alpint < alphaMax ? yldred*pow(1.+beta*a->alpint,npow) : yldredMin ;
}

// Get derivative of sqrt(2./3.)*yield with respect to lambda for plane strain and 3D
// ... and using dep/dlambda = sqrt(2./3.)
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)lamda or depdot/dlambda = sqrt(2./3.)/delTime
double NonlinearHardening::GetKPrime(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	return a->alpint < alphaMax ? TWOTHIRDS*yldred*beta*npow*pow(1.+beta*a->alpint,npow-1) : 0. ;
}

// Get derivative of (1./3.)*yield^2 with respect to lambda for plane stress only
// ... and using dep/dlambda = sqrt(2./3.)*fnp1 where ep=alpint
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)*lambda*fnp1 or depdot/dlambda = sqrt(2./3.)*fnp1/delTime
// Also equal to sqrt(2./3.)*GetYield()*GetKPrime()*fnp1, but in separate call for efficiency
double NonlinearHardening::GetK2Prime(MPMBase *mptr,double fnp1,double delTime,HardeningAlpha *a,void *properties) const
{
	return a->alpint < alphaMax ? SQRT_EIGHT27THS*yldred*yldred*beta*npow*pow(1.+beta*a->alpint,2.*npow-1)*fnp1 : 0. ;
}

#pragma mark NonlinearHardening::Accessors

// hardening law name
const char *NonlinearHardening::GetHardeningLawName(void) const { return "Nonlinear hardening"; }

