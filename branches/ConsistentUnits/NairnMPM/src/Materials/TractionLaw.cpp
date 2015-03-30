/********************************************************************************
    TractionLaw.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Feb 22 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/TractionLaw.hpp"
#include "Cracks/CrackSegment.hpp"
#include "System/ArchiveData.hpp"
#include "System/UnitsController.hpp"

#pragma mark TractionLaw::Constructors and Destructors

// Constructors with arguments 
TractionLaw::TractionLaw(char *matName) : MaterialBase(matName)
{
	stress1=-1;			// sigmaI
	stress2=-1;			// sigmaII
}

#pragma mark TractionLaw::Initialization

// do not print the base class transport properties
void TractionLaw::PrintTransportProperties(void) const {}

// read peak stress in all traction laws
char *TractionLaw::InputMaterialProperty(char *xName,int &input,double &gScaling)
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

// report debond in same format for all cohesive laws
// dtime in sec, cs if the debonded segment, fractionI is fraction mode I at time of debond
void TractionLaw::ReportDebond(double dtime,CrackSegment *cs,double fractionI,double Gtotal)
{
	archiver->IncrementPropagationCounter();
	cout << "# Debond: t=" << dtime*UnitsController::Scaling(1000.) << " (x,y) = (" << cs->x << "," << cs->y << ")"
			<< " GI(%) = " << 100.*fractionI << " G = "
			<< Gtotal*UnitsController::Scaling(0.001) << endl;
}

#pragma mark TractionLaw::Traction Law

// Traction law - find traction force as traction pressure*area
void TractionLaw::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,double dx,double dy,double area)
{	cs->tract.x=0.;
	cs->tract.y=0.;
}

// Find recoverable energy in the traction law (Subclass must override)
// Units are work/area or force/length
double TractionLaw::CrackTractionEnergy(CrackSegment *cs,double nCod,double tCod,bool fullEnergy)
{	return 0.;
}

#pragma mark TractionLaw::Accessors

// required accessors
double TractionLaw::WaveSpeed(bool threeD,MPMBase *mptr) const { return 1.e-12; }

// return material type
const char *TractionLaw::MaterialType(void) const { return "Crack Traction Law"; }

// check if traction law material
bool TractionLaw::isTractionLaw(void) const { return TRUE; }

