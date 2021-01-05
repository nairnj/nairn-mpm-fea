/********************************************************************************
    ExponentialTraction.cpp
    nairn-mpm-fea

    Created by John Nairn on 12/22/2020.
    Copyright (c) 2020 John A. Nairn, All rights reserved.
 
    This material used the scaled exponential law allowed in Abaqus
********************************************************************************/

#include "stdafx.h"
#include "Materials/ExponentialTraction.hpp"

extern double mtime;

#pragma mark ExponentialTraction::Constructors and Destructors

// Constructor
ExponentialTraction::ExponentialTraction(char *matName,int matID) : CohesiveZone(matName,matID)
{
    // user must provide positive values
    alphaI = -1.;
    alphaII = -1.;
}

#pragma mark ExponentialTraction::Initialization

// Read properteis (read read in super classes)
char *ExponentialTraction::InputTractionLawProperty(char *xName,int &input,double &gScaling)
{
    if(strcmp(xName,"alphaI")==0)
    {   input=DOUBLE_NUM;
        return((char *)&alphaI);
    }
    
    else if(strcmp(xName,"alphaII")==0)
    {   input=DOUBLE_NUM;
        return((char *)&alphaII);
    }
    
    return CohesiveZone::InputTractionLawProperty(xName,input,gScaling);
}

// Calculate properties used in analyses - here triangular law
// Do mode I and mode II separately
const char *ExponentialTraction::VerifyAndLoadProperties(int np)
{
    if(alphaI<=0.)
        return "Exponential cohesive laws must specify alphaI > 0";
    expalphaI = exp(-alphaI);
    romexpalphaI = 1./(1.-expalphaI);
    double fealphaI = 2.*(1./alphaI - expalphaI*romexpalphaI);
    const char *msg=SetExponentialTractionLaw(stress1,kI1,delIc,JIc,umidI,fealphaI);
    if(msg!=NULL) return msg;
    //cout << "# I: " << expalphaI << "," << romexpalphaI << "," << fealphaI << endl;
    
    if(alphaII<=0.)
        return "Exponential cohesive laws must specify alphaII > 0";
    expalphaII = exp(-alphaII);
    romexpalphaII = 1./(1.-expalphaII);
    double fealphaII = 2.*(1./alphaII - expalphaII*romexpalphaII);
    msg=SetExponentialTractionLaw(stress2,kII1,delIIc,JIIc,umidII,fealphaII);
    if(msg!=NULL) return msg;
    //cout << "# I: " << expalphaI << "," << romexpalphaI << "," << fealphaI << endl;

    // go to parent
    return TractionLaw::VerifyAndLoadProperties(np);
}

// print to output window
void ExponentialTraction::PrintMechanicalProperties(void) const
{
    PrintSawToothModel("I",JIc,stress1,delIc,kI1,umidI,alphaI);
    PrintSawToothModel("II",JIIc,stress2,delIIc,kII1,umidII,alphaII);
    PrintProperty("n",nmix,"");
    cout << endl;
}

#pragma mark CohesiveZone::Basic Functions

// Return the strength for mode and current delta
// delta must be in valid range
double ExponentialTraction::Strength(int mode,double delta)
{   double arg;
    if(mode==1)
    {   arg = exp(-alphaI*(delta-umidI)/(delIc-umidI));    // = 1 when still elastic
        return stress1*romexpalphaI*(arg-expalphaI);
    }
    arg = exp(-alphaII*(delta-umidII)/(delIIc-umidII));    // = 1 when still elastic
    return stress2*romexpalphaII*(arg-expalphaII);
}

// Return area under the coshesive law up to u (only used in J integral)
// Assumes ue <= u <= uc
double ExponentialTraction::WorkEnergy(int mode,double u)
{   double arg;
    if(mode==1)
    {   arg = exp(-alphaI*(u-umidI)/(delIc-umidI));
        return 0.5*stress1*( umidI + 2.*romexpalphaI*( (delIc-umidI)*(1.-arg)/alphaI
                                   - (u-umidI)*expalphaI ) );
    }
    arg = exp(-alphaII*(u-umidII)/(delIIc-umidII));
    return 0.5*stress2*( umidII + 2.*romexpalphaII*( (delIIc-umidII)*(1.-arg)/alphaII
                                    - (u-umidII)*expalphaII ) );
}

// Return dissipated energy up to delta.
// delta must be in valid range (delta>umid)
double ExponentialTraction::DissipatedEnergy(int mode,double delta)
{   double arg;
    if(mode==1)
    {   arg = exp(-alphaI*(delta-umidI)/(delIc-umidI));
        return 0.5*stress1*( umidI + romexpalphaI*( 2.*(delIc-umidI)*(1.-arg)/alphaI
                    - delta*(expalphaI+arg) + 2.*umidI*expalphaI ) );
    }
    arg = exp(-alphaII*(delta-umidII)/(delIIc-umidII));
    return 0.5*stress2*( umidII + romexpalphaII*( 2.*(delIIc-umidII)*(1.-arg)/alphaII
                    - delta*(expalphaII+arg) + 2.*umidII*expalphaII ) );
}

// Get D from delta (needed for MixedModeTraction)
// delta must be in valid range
double ExponentialTraction::GetDFromDelta(int mode,double delta)
{   double arg;
    if(mode==1)
    {   arg = exp(-alphaI*(delta-umidI)/(delIc-umidI));
        return 1. - (umidI/delta)*(1.-romexpalphaI*(1.-arg));
    }
    arg = exp(-alphaII*(delta-umidII)/(delIIc-umidII));
    return 1. - (umidII/delta)*(1.-romexpalphaII*(1.-arg));
}

// Get delta from D (needed for MixedModeTraction)
double ExponentialTraction::GetDeltaFromD(int mode,double D)
{   double uea;
    if(mode==1)
    {   uea = umidI*romexpalphaI*expalphaI/(1.-D);
        return (delIc-umidI)*gsl_sf_lambert_W0(alphaI*uea*exp(alphaI*(delIc+uea)/(delIc-umidI))/(delIc-umidI))/alphaI - uea;
    }
    uea = umidII*romexpalphaII*expalphaII/(1.-D);
    return (delIIc-umidII)*gsl_sf_lambert_W0(alphaII*uea*exp(alphaII*(delIIc+uea)/(delIIc-umidII))/(delIIc-umidII))/alphaII - uea;
}

// Return dissipation rate function (phi(delta))
// delta must be in valid range
double ExponentialTraction::DissipationRate(int mode,double delta)
{   double arg;
    if(mode==1)
    {   arg = exp(-alphaI*(delta-umidI)/(delIc-umidI));    // = 1 when still elastic
        return stress1*romexpalphaI*(arg*(1+alphaI*delta/(delIc-umidI))-expalphaI);
    }
    arg = exp(-alphaII*(delta-umidII)/(delIIc-umidII));    // = 1 when still elastic
    return stress2*romexpalphaII*(arg*(1+alphaII*delta/(delIIc-umidII))-expalphaII);
}

// Return uc*dD/ddelta or ratio of damage parameter evolution to delta evolution
// delta must be in valid range, output is dimensionless
// Only needed if MixedModeTraction calls it
double ExponentialTraction::RatioFunction(int mode,double delta)
{   double varphi = ExponentialTraction::DissipationRate(mode,delta);
    return mode==1 ? varphi*umidI/(kI1*delta*delta) : varphi*umidII/(kII1*delta*delta) ;
}

#pragma mark ExponentialTraction::Accessors

// return material type
const char *ExponentialTraction::MaterialType(void) const { return "Exponential Cohesive Law"; }

