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

// Constructors
HillPlastic::HillPlastic()
{
}

// Constructors
HillPlastic::HillPlastic(char *matName) : AnisoPlasticity(matName)
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

// hardening history - equivalent plastic strain (absolute strain)
double HillPlastic::GetHistory(int num,char *historyPtr) const
{
    double history=0.;
	if(num==1)
	{	double *h=(double *)historyPtr;
		history=*h;
	}
    return history;
}

#pragma mark HillPlastic:Hardening Terms

// Load current internal variables into local alpha variables
void HillPlastic::UpdateTrialAlpha(MPMBase *mptr,int np,AnisoHardProperties *p,int offset) const
{	p->aint = mptr->GetHistoryDble(0,offset);
}

// Update alpha: Here dalpha = alpha0 - lambdak h = alpha0 + lambdak*(-h)
void HillPlastic::UpdateTrialAlpha(MPMBase *mptr,int np,double lambdak,AnisoHardProperties *p,int offset) const
{	p->aint = mptr->GetHistoryDble(0,offset) + lambdak*p->minush;
}

// Return yield stress for current conditions (p->aint for cum. plastic strain and dalpha/delTime for plastic strain rate)
double HillPlastic::GetYield(AnisoHardProperties *p) const
{
	return 1. + Khard*pow(p->aint,nhard);
}

// Hardening term : find df^(alpha) . h = (-g'(alpha))*(-h)
// Here g(alpha) = 1 + Khard alpha^nhard
//    and therefore g'(alpha) = nhard Khard alpha^(nhard-1) or just Khard if nhard=1
double HillPlastic::GetDfAlphaDotH(MPMBase *mptr,int np,AnisoHardProperties *p) const
{	return DbleEqual(nhard,1.) ? Khard*p->minush :
	Khard*nhard*pow(p->aint,nhard-1)*p->minush ;
}

// transfer final alpha variables to the material point
void HillPlastic::UpdatePlasticInternal(MPMBase *mptr,int np,AnisoHardProperties *p,int offset) const
{	mptr->SetHistoryDble(0,p->aint,offset);
}


#pragma mark HillPlastic:Accessors

// return material type
const char *HillPlastic::MaterialType(void) const { return "Elastic-Plastic Hill Material"; }

