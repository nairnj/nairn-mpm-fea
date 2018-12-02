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
        return((char *)&nhard);
    }
	
    return(AnisoPlasticity::InputMaterialProperty(xName,input,gScaling));
}

// print to output window
void HillPlastic::PrintYieldProperties(void) const
{	
    AnisoPlasticity::PrintYieldProperties();
	PrintProperty("K",Khard,"");
	PrintProperty("n",nhard,"");
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
{	return 1. + Khard*pow(p->aint,nhard);
}

// Get g'(alpha)
double HillPlastic::GetGPrime(AnisoPlasticProperties *p) const
{	return DbleEqual(nhard,1.) ? Khard : Khard*nhard*pow(p->aint,nhard-1) ;
}


#pragma mark HillPlastic:Accessors

// return material type
const char *HillPlastic::MaterialType(void) const { return "Elastic-Plastic Hill Material"; }

