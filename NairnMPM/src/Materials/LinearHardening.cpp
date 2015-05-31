/********************************************************************************
    LinearHardening.hpp
    nairn-mpm-fea

    Created by John Nairn, 1/17/2103
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Hardening Law is
        sigma = yield (1 + K ep)
        Use zero to get elastic plastic
 
    Previous code used yield + Ep ep   or   yield*K = Ep
        It also allowed entry of ET which was coverted to Ep=E ET/(E-ET)
        or 1/K = yield((1/ET) - (1/E)). The switch had to be made because
        this law can be used with materials that do not define E
********************************************************************************/

#include "Materials/LinearHardening.hpp"
#include "System/UnitsController.hpp"

#pragma mark LinearHardening::Constructors and Destructors

// Constructors
LinearHardening::LinearHardening() {}

// Constructors
LinearHardening::LinearHardening(MaterialBase *pair) : HardeningLawBase(pair)
{
    beta = 0.;
    Ep = -1.;
}

#pragma mark LinearHardening::Initialize

// Read hardening law properties
char *LinearHardening::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    // dimensionless coefficient for hardening
    if(strcmp(xName,"Khard")==0)
    {   input=DOUBLE_NUM;
        return (char *)&beta;
    }

    else if(strcmp(xName,"Ep")==0)
    {   input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&Ep,gScaling,1.e6);
    }
    
    return HardeningLawBase::InputMaterialProperty(xName,input,gScaling);
}

// get reduced stress than done
const char *LinearHardening::VerifyAndLoadProperties(int np)
{
    // call first to get reduced yield stress
    HardeningLawBase::VerifyAndLoadProperties(np);
    
    // Use Ep if it was entered and is non-negative
    if(Ep >= 0.)
		beta = Ep/yield;
	else
		Ep = beta*yield;
    
	// save some multiplies (it is reduced plastic modulus)
    Epred = yldred*beta;
	
	// base call above never has an error
	return NULL;
}

// print just yield properties to output window
void LinearHardening::PrintYieldProperties(void) const
{
    cout << GetHardeningLawName() << endl;
    MaterialBase::PrintProperty("yld",yield*UnitsController::Scaling(1.e-6),"");
    MaterialBase::PrintProperty("K",beta,"");
    MaterialBase::PrintProperty("Ep",Ep*UnitsController::Scaling(1.e-6),"");
    cout << endl;
}

#pragma mark LinearHardening::Law Methods

// Return yield stress for current conditions (alpint for cum. plastic strain and dalpha/delTime for plastic strain rate)
double LinearHardening::GetYield(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	return yldred + Epred*a->alpint;
}

// Return (K(alpha)-K(0)), which is used in dissipated energy calculation
double LinearHardening::GetYieldIncrement(MPMBase *mptr,int np,double delTime,HardeningAlpha *a,void *properties) const
{
	return Epred*a->alpint;
}

// Get derivative of sqrt(2./3.)*yield with respect to lambda for plane strain and 3D
// ... and using dep/dlambda = sqrt(2./3.)
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)lamda or depdot/dlambda = sqrt(2./3.)/delTime
double LinearHardening::GetKPrime(MPMBase *mptr,int np,double delTime,HardeningAlpha *,void *properties) const
{
	return TWOTHIRDS*Epred;
}

// Get derivative of (1./3.)*yield^2 with respect to lambda for plane stress only
// ... and using dep/dlambda = sqrt(2./3.)*fnp1 where ep=alpint
// ... and epdot = dalpha/delTime with dalpha = sqrt(2./3.)*lambda*fnp1 or depdot/dlambda = sqrt(2./3.)*fnp1/delTime
// Also equal to sqrt(2./3.)*GetYield()*GetKPrime()*fnp1, but in separate call for efficiency
double LinearHardening::GetK2Prime(MPMBase *mptr,double fnp1,double delTime,HardeningAlpha *a,void *properties) const
{
	return SQRT_EIGHT27THS*(yldred + Epred*a->alpint)*Epred*fnp1;
}

#pragma mark HardeningLawBase::Return Mapping

// Linear law can do return mapping analytically, except for plane stress
double LinearHardening::SolveForLambdaBracketed(MPMBase *mptr,int np,double strial,Tensor *stk,double Gred,
												double psKred,double Ptrial,double delTime,HardeningAlpha *a,void *p) const
{
	// plane stress is numerical
	if(np==PLANE_STRESS_MPM)
    {   // The unbracketed one is faster and seems stable for this hardening law
        return HardeningLawBase::SolveForLambda(mptr,np,strial,stk,Gred,psKred,Ptrial,delTime,a,p);
    }
    
	// closed form for plane strain and 3D
	double lambdak = (strial - SQRT_TWOTHIRDS*(yldred + Epred*a->alpint))/(2.*(Gred + Epred/3.));
	UpdateTrialAlpha(mptr,np,lambdak,(double)1.,a);
	return lambdak;
}

#pragma mark LinearHardening::Accessors

// hardening law name
const char *LinearHardening::GetHardeningLawName(void) const { return "Linear hardening"; }

