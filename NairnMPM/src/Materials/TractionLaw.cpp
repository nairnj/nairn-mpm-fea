/********************************************************************************
    TractionLaw.cpp
    NairnMPM
    
    Created by John Nairn on Feb 22 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/TractionLaw.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include "Cracks/CrackSegment.hpp"

#pragma mark TractionLaw::Constructors and Destructors

// Constructors with arguments 
TractionLaw::TractionLaw(char *matName) : MaterialBase(matName)
{
	stress1=-1;			// sigmaI
	stress2=-1;			// sigmaII
}

#pragma mark TractionLaw::Initialization

// do not print the base class transport properties
void TractionLaw::PrintTransportProperties(void) {}

// read peak stress in all traction laws
char *TractionLaw::InputMat(char *xName,int &input)
{
    if(strcmp(xName,"sigmaI")==0)
	{	input=DOUBLE_NUM;
        return((char *)&stress1);
    }

	else if(strcmp(xName,"sigmaII")==0)
	{	input=DOUBLE_NUM;
		return((char *)&stress2);
	}
	
    return MaterialBase::InputMat(xName,input);
}

// do not need to call base material class methods
const char *TractionLaw::VerifyProperties(int np) { return NULL; }

#pragma mark TractionLaw::Traction Law

// Traction law - subclass must override
void TractionLaw::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,double dx,double dy,double area)
{	cs->tract.x=0.;
	cs->tract.y=0.;
}

// return recoverable energy (I think this is correct) in the traction law
// units of N/mm. Subclass must override
double TractionLaw::CrackTractionEnergy(CrackSegment *cs,double nCod,double tCod,bool fullEnergy)
{	return 0.;
}

// define empty methods for constuitive laws that are not used
void TractionLaw::MPMConstLaw(MPMBase *,double,double,double,double,double,int) {}
void TractionLaw::MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int) {}

#pragma mark TractionLaw::Accessors

// required accessors
double TractionLaw::WaveSpeed(bool threeD) { return 1.e-12; }

// return material type
const char *TractionLaw::MaterialType(void) { return "Crack Traction Law"; }

// check if traciton law material
bool TractionLaw::isTractionLaw(void) { return TRUE; }

