/********************************************************************************
    DugdaleMat.cpp
    NairnMPM
    
    Created by John Nairn on Wed Sep 04 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/DugdaleMat.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"

#pragma mark DugdaleMat::Constructors and Destructors

// Constructors
DugdaleMat::DugdaleMat() {}

// Constructors
DugdaleMat::DugdaleMat(char *matName) : VonMisesHardening(matName)
{
    matType=DUGDALE;
}

#pragma mark DugdaleMat::Methods

// For Dugdale linear hardening material, plastic potential 
// f=s_y(material axes)^2-yield^2-2 Ep Wp
// Ep - plastic modulus; Wp - plastic work.
double DugdaleMat::GetF(MPMBase *mptr,Tensor *st,int np)
{
	double angle=mptr->GetRotation();
    double Wp=mptr->GetPlastEnergy();
	double c=cos(angle);
	double c2=c*c;
	double s=sin(angle);
	double s2=s*s;
	// clockwise rotation from analysis to material axes
	double syrot=st->xx*s2 + st->yy*c2 + 2.*st->xy*c*s;
	return syrot*syrot - yldred*yldred - 2.*Epred*Wp;
}

// return derivatives of the yield function wrt all stresses
void DugdaleMat::GetDfDsigma(MPMBase *mptr,Tensor *st0,int np)
{
	double angle=mptr->GetRotation();
	double c=cos(angle);
	double c2=c*c;
	double s=sin(angle);
	double s2=s*s;
	// clockwise rotation from analysis to material axes
	double syrot=(st0->xx*s2 + st0->yy*c2 + 2.*st0->xy*c*s);
	dfdsxx=2.*s2*syrot;
	dfdsyy=2.*c2*syrot;
	dfdtxy=4.*c*s*syrot;
	dfdszz=0.;
}

#pragma mark DugdaleMat:Accessors

// return material type
char *DugdaleMat::MaterialType(void) { return "Isotropic Dugdale"; }
