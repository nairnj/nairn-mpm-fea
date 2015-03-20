/********************************************************************************
    NonlinearHardening.hpp
    nairn-mpm-fea

    Created by John Nairn, 1/17/2103
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Hardening Law is
        sigma = yield (1 + beta ep)^npow
********************************************************************************/

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
}

#pragma mark LinearHardening::Initialize

// Read hardening law properties
char *NonlinearHardening::InputMaterialProperty(char *xName,int &input,double &gScaling)
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
    
    return HardeningLawBase::InputMaterialProperty(xName,input,gScaling);
}

// print just yield properties to output window
void NonlinearHardening::PrintYieldProperties(void) const
{
    cout << GetHardeningLawName() << endl;
    MaterialBase::PrintProperty("yld",yield*UnitsController::Scaling(1.e-6),"");
    MaterialBase::PrintProperty("beta",beta,"");
    MaterialBase::PrintProperty("n",npow,"");
    cout << endl;
}

#pragma mark NonlinearHardening::Law Methods

// Return yield stress for current conditions (alpint for cum. plastic strain and dalpha/delTime for plastic strain rate)
double NonlinearHardening::GetYield(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{   
	return yldred*pow(1.+beta*a->alpint,npow) ;
}

// Get derivative of sqrt(2./3.)*yield with respect to lambda for plane strain and 3D
// ... and using dep/dlambda = sqrt(2./3.)
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)lamda or depdot/dlambda = sqrt(2./3.)/delTime
double NonlinearHardening::GetKPrime(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	return TWOTHIRDS*yldred*beta*npow*pow(1.+beta*a->alpint,npow-1) ;
}

// Get derivative of (1./3.)*yield^2 with respect to lambda for plane stress only
// ... and using dep/dlambda = sqrt(2./3.)*fnp1 where ep=alpint
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)*lambda*fnp1 or depdot/dlambda = sqrt(2./3.)*fnp1/delTime
// Also equal to sqrt(2./3.)*GetYield()*GetKPrime()*fnp1, but in separate call for efficiency
double NonlinearHardening::GetK2Prime(MPMBase *mptr,double fnp1,double delTime,HardeningAlpha *a,void *properties) const
{
	return SQRT_EIGHT27THS*yldred*yldred*beta*npow*pow(1.+beta*a->alpint,2.*npow-1)*fnp1;
}

#pragma mark NonlinearHardening::Accessors

// hardening law name
const char *NonlinearHardening::GetHardeningLawName(void) const { return "Nonlinear hardening"; }

