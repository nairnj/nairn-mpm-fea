/********************************************************************************
    HillPlastic.cpp
    nairn-mpm-fea
    
    Created by John Nairn, June 12, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		AnisoPlasticity.hpp (Orthotropic.hpp, TranIsotropic.hpp, MaterialBase.hpp)
********************************************************************************/

#include "stdafx.h"
#include "HillPlastic.hpp"
#include "MPM_Classes/MPMBase.hpp"

#pragma mark HillPlastic::Constructors and Destructors

// Constructor
HillPlastic::HillPlastic(char *matName,int matID) : AnisoPlasticity(matName,matID)
{
    // default value of hardening (elastic plastic)
	Khard=0.;
	nhard=1.;
    exphard=0.;
    hardStyle = AP_LINEAR;
}

#pragma mark HillPlastic::Initialize

// Read material properties
char *HillPlastic::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    // Khard: hardening parameter
    if(strcmp(xName,"Khard")==0)
    {   input=DOUBLE_NUM;
        return((char *)&Khard);
    }

    // nhard - hardening exponent
    else if(strcmp(xName,"nhard")==0)
    {   input=DOUBLE_NUM;
        hardStyle = AP_UNKNOWN;
        return((char *)&nhard);
    }
	
    // exphard - hardening exponential rate
    else if(strcmp(xName,"exphard")==0)
    {   input=DOUBLE_NUM;
        hardStyle = AP_EXPONENTIAL;
        return((char *)&exphard);
    }

    return(AnisoPlasticity::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
// NOTE: This code duplicated in OrthoPlasticSoftening. Keep them in sync
const char *HillPlastic::VerifyAndLoadProperties(int np)
{
    if(hardStyle==AP_UNKNOWN)
    {   hardStyle = nhard<=0. ? AP_NONLINEAR1 : AP_NONLINEAR2;
        nhard = fabs(nhard);
    }
    else if(hardStyle==AP_EXPONENTIAL)
    {   if(exphard==0.)
            return "The exphard parameter cannot be zero.";
        Kexp = Khard/exphard;
    }
    
    // call super class
    return AnisoPlasticity::VerifyAndLoadProperties(np);
}

// print to output window
void HillPlastic::PrintYieldProperties(void) const
{	
    AnisoPlasticity::PrintAPYieldProperties(syxx,syyy,syzz,tyyz,tyxz,tyxy);
    cout << "Relative yield = ";
    switch(hardStyle)
    {   case AP_LINEAR:
            cout << "1 + " << Khard << "*alpha";
            break;
        case AP_NONLINEAR1:
            cout << "(1 + " << Khard << "*alpha)^" << nhard;
            break;
        case AP_EXPONENTIAL:
            if(exphard>0.)
                cout << "1 + (" << Kexp << ")*(1-exp[-" << exphard << "*alpha])";
            else
                cout << "1 + (" << Kexp << ")*(1-exp[" << fabs(exphard) << "*alpha])";
            break;
        default:
            cout << "1 + " << Khard << "*alpha^" << nhard;
            break;
    }
    cout << endl;
}

#pragma mark HillPlastic:History Data Methods

// history is cumulative strain
char *HillPlastic::InitHistoryData(char *pchr,MPMBase *mptr)
{
	double *p = CreateAndZeroDoubles(pchr,1);
	return (char *)p;
}

// Number of history variables
int HillPlastic::NumberOfHistoryDoubles(void) const { return 1; }

#pragma mark HillPlastic:Hardening Terms

// Return yield stress for current conditions (p->aint for cum. plastic strain and dalpha/delTime for plastic strain rate)
double HillPlastic::GetYield(AnisoPlasticProperties *p) const
{   switch(hardStyle)
    {   case AP_LINEAR:
            return 1. + Khard*p->aint;
        case AP_NONLINEAR1:
            return pow(1. + Khard*p->aint,nhard);
        case AP_EXPONENTIAL:
            return 1. + Kexp*(1.-exp(-exphard*p->aint));
        default:
            return 1. + Khard*pow(p->aint,nhard);
     }
}

// Get g'(alpha)
double HillPlastic::GetGPrime(AnisoPlasticProperties *p) const
{	switch(hardStyle)
    {   case AP_LINEAR:
            return Khard;
        case AP_NONLINEAR1:
            return Khard*nhard*pow(1. + Khard*p->aint,nhard-1.);
        case AP_EXPONENTIAL:
            return Khard*exp(-exphard*p->aint);
        default:
            if(nhard<1. && p->aint<ALPHA_EPS)
                return Khard*nhard*pow(ALPHA_EPS,nhard-1.);
            else
                return Khard*nhard*pow(p->aint,nhard-1.);
    }
}

#pragma mark HillPlastic:Accessors

// return material type
const char *HillPlastic::MaterialType(void) const { return "Elastic-Plastic Hill Material"; }

