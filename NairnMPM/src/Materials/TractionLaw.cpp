/********************************************************************************
    TractionLaw.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Feb 22 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/TractionLaw.hpp"
#include "Cracks/CrackSegment.hpp"
#include "System/ArchiveData.hpp"
#include "System/UnitsController.hpp"

#pragma mark TractionLaw::Constructors and Destructors

// Constructor 
TractionLaw::TractionLaw(char *matName,int matID) : MaterialBase(matName,matID)
{
	stress1=-1;			// sigmaI
	stress2=-1;			// sigmaII
    
    numTractionHistory = 0; // set to number of history doubles used by the law
}

#pragma mark TractionLaw::Initialization

// do not print the base class transport properties
void TractionLaw::PrintTransportProperties(void) const {}

// read peak stress in all traction laws
char *TractionLaw::InputMaterialProperty(char *xName,int &input,double &gScaling)
{	return InputTractionLawProperty(xName,input,gScaling);
}
char *TractionLaw::InputTractionLawProperty(char *xName,int &input,double &gScaling)
{
    if(strcmp(xName,"sigmaI")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&stress1,gScaling,1.e6);
    }

	else if(strcmp(xName,"sigmaII")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&stress2,gScaling,1.e6);
	}
	
    return MaterialBase::InputMaterialProperty(xName,input,gScaling);
}

// do not need to call base material class methods
const char *TractionLaw::VerifyAndLoadProperties(int np) { return NULL; }

// pack label for printing (label must be long enough)
char *TractionLaw::PackLabel(char *label,const char *prefix,const char *mode,const char *suffix) const
{   strcpy(label,prefix);
    strcat(label,mode);
    strcat(label,suffix);
    return label;
}

// report debond in same format for all cohesive laws
// dtime in sec, cs if the debonded segment, fractionI is fraction mode I at time of debond
void TractionLaw::ReportDebond(double dtime,CrackSegment *cs,double fractionI,double Gtotal)
{
	archiver->IncrementPropagationCounter();
	cout << "# Debond: t=" << dtime*UnitsController::Scaling(1000.) << " (x,y) = (" << cs->cp.x << "," << cs->cp.y << ")"
			<< " GI(%) = " << 100.*fractionI << " G = "
			<< Gtotal*UnitsController::Scaling(0.001) << endl;
}

// evaluate pressure at current time
// (Don't call in parallel code due to function)
void TractionLaw::CalculateTimeFunction(void) {}

#pragma mark TractionLaw::Traction Law

// Traction law - find traction force as traction pressure*area
void TractionLaw::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,Vector *n,Vector *t,double area)
{	cs->tract.x=0.;
	cs->tract.y=0.;
	cs->tract.z=0.;
}

// Return current traction law strain energy (Int T.du).
//	This energy is needed for J integral (and only used in J Integral)
// units of F/L
double TractionLaw::CrackWorkEnergy(CrackSegment *cs,double nCod,double tCod)
{	return 0.;
}

// Return mode I and II energy that has been released by current segment.
void TractionLaw::CrackDissipatedEnergy(CrackSegment *cs,double &GI,double &GII)
{	GI=0.;
	GII=0.;
}

#pragma mark TractionLaw::Basic Functions

// Warning - this was added after mch coding was done. Not all
// traction laws use them

// Return the strength for mode and current delta
// delta must be in valid range
double TractionLaw::Strength(int mode,double delta) { return 0.; }

// Return area under the coshesive law up to u (<= uc)
// Only used in J integral
double TractionLaw::WorkEnergy(int mode,double u) { return 0.; }

// Return dissipated energy up to delta
// delta must be in valid range
double TractionLaw::DissipatedEnergy(int mode,double delta) { return 0.; }

// Return the derivative of strength with repect to delta
// delta must be in valid range
double TractionLaw::StrengthPrime(int mode,double delta) { return 0.; }

// Return dissipaation rate function (phi(delta))
// delta must be in valid range
// Only needed if MixedModeTraction calls it
double TractionLaw::DissipationRate(int mode,double delta) { return 0.; }

// Return uc*dD/ddelta or ratio of damage parameter evolution to delta evolution
// delta must be in valid range, output is dimensionless
// Only needed if MixedModeTraction calls it
double TractionLaw::RatioFunction(int mode,double delta) { return 0.; }

// Get D from delta (needed for MixedModeTraction)
// delta must be in valid range
double TractionLaw::GetDFromDelta(int mode,double delta) { return 1.e10; }

// Get delta from D (needed for MixedModeTraction)
double TractionLaw::GetDeltaFromD(int mode,double D) { return 1.; }

#pragma mark TractionLaw::Cubic Law Specific

// set traction law on initialization
//        smax is peak stess (F/L^2)
//        umax is failure displacement (mm)
//        k1 is slope (an output, not an input) (F/L^5)
//        G is toughness (F/L) = (9/16) umax smax
const char *TractionLaw::SetCubicTractionLaw(double &smax,double &umax,double &G,double &k1)
{
    // specify smax and G
    if(umax<0.)
    {    if(smax<0. || G<0.)
            return "Too few cubic traction law properties were supplied.";
        umax = 16.*G/(9.*smax);
    }
    
    // specify umax and G
    else if(smax<0.)
    {    if(G<0.)
            return "Too few cubic traction law properties were supplied.";
        smax = 16.*G/(9.*umax);
    }
    
    // specify umax and smax
    else if(G<0.)
    {    G = 9.*umax*smax/16.;
    }
    
    // specified them all, which is an error
    else
    {    return "Must supply exactly two of delIc, sigmaI, JIc and exactly two of delIIc, sigmaII, JIIc.";
    }
    
    // stress prefactor to get force
    k1 = 27.*smax/(4.*umax*umax*umax);
    
    return NULL;
}

// print one cubic mode law
void TractionLaw::PrintCubicModel(const char *mode,double Jc,double sigc,double uc,double ke) const
{
    char label[10];
    PrintProperty(PackLabel(label,"G",mode,"c"),Jc*UnitsController::Scaling(0.001),UnitsController::Label(ERR_UNITS));
    PrintProperty(PackLabel(label,"sig",mode,""),sigc*UnitsController::Scaling(1.e-6),UnitsController::Label(PRESSURE_UNITS));
    PrintProperty(PackLabel(label,"u",mode,""),uc,UnitsController::Label(CULENGTH_UNITS));
    PrintProperty(PackLabel(label,"k",mode,""),ke*delIc*delIc*UnitsController::Scaling(1.e-6),UnitsController::Label(TRACTIONSLOPE_UNITS));
    cout <<  endl;
}

#pragma mark TractionLaw::Accessors

// required accessors
double TractionLaw::WaveSpeed(bool threeD,MPMBase *mptr) const { return 1.e-12; }

// return material type
const char *TractionLaw::MaterialType(void) const { return "Crack Traction Law"; }

// check if traction law material
int TractionLaw::MaterialStyle(void) const { return TRACTION_MAT; }

// Number of history variables (used by some materials when creating array of doubles)
int TractionLaw::NumberOfHistoryDoubles(void) const { return numTractionHistory; }


