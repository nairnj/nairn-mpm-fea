/********************************************************************************
	TrilinearTraction.cpp
	nairn-mpm-fea

	Created by John Nairn on 9/17/10.
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	The traction law is trilinear with break points at (u1,s1)
	and (u2,s2). u1 may be zero for bilinear softening.

	s1|     +
	  |    + +
	  |    +  +
	  |   +    +
	  |   +     +
	  |  +       +
	  |  +        +
	s2| +           +++
	  | +              +++
	  |±__________________±±±______
	         u1     u2         u3
 
	Mode I: s1=stress1, u1=umidI, s2=sI2, u2=uI2, u3=delIc, area=JIc
	Mode II: s1=stress2, u1=umidII, s2=sII2, u2=uII2, u3=delIIc, area=JIc
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Materials/TrilinearTraction.hpp"
#include "Cracks/CrackSegment.hpp"
#include "System/UnitsController.hpp"

extern double mtime;

#pragma mark TrilinearTraction::Constructors and Destructors

// Constructor 
TrilinearTraction::TrilinearTraction(char *matName,int matID) : ExponentialTraction(matName,matID)
{
	// mode I cohesive law (all others set to -1 in superclasses)
	// others are: stress1,delIc,JIc,kI1,umid1
	sI2=-1.;				// second stress
	uI2=-1.;				// second displacement point
	
	// mode II cohesive law (all others set to -1 in superclasses)
	// others are: stress2,delIIc,JIIc,kII1,umidII
	sII2=-1.;				// second stress
	uII2=-1.;				// second displacement point
}

#pragma mark TrilinearTraction::Initialization

// no properties to read
char *TrilinearTraction::InputTractionLawProperty(char *xName,int &input,double &gScaling)
{
	if(strcmp(xName,"sigmaI2")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&sI2,gScaling,1.e6);
	}
	
	else if(strcmp(xName,"sigmaII2")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&sII2,gScaling,1.e6);
	}
	
	else if(strcmp(xName,"delpkI2")==0)
	{	input=DOUBLE_NUM;
		return((char *)&uI2);
	}
	
	else if(strcmp(xName,"delpkII2")==0)
	{	input=DOUBLE_NUM;
		return((char *)&uII2);
	}
	
    return ExponentialTraction::InputTractionLawProperty(xName,input,gScaling);
}

// Calculate properties used in analyses - here trilinear law
const char *TrilinearTraction::VerifyAndLoadProperties(int np)
{
	const char *msg=SetTLTractionLaw(stress1,kI1,umidI,sI2,uI2,delIc,JIc,break1is2I,JI_1c,JI_2c);
	if(msg!=NULL) return msg;
	
	msg=SetTLTractionLaw(stress2,kII1,umidII,sII2,uII2,delIIc,JIIc,break1is2II,JII_1c,JII_2c);
	if(msg!=NULL) return msg;
	
	return TractionLaw::VerifyAndLoadProperties(np);
}

// print to output window
void TrilinearTraction::PrintMechanicalProperties(void) const
{
    PrintTriLinearModel("I",JIc,stress1,kI1,umidI,sI2,uI2,delIc);
    PrintTriLinearModel("II",JIIc,stress2,kII1,umidII,sII2,uII2,delIIc);
}

#pragma mark CohesiveZone::Basic Functions

// Return the strength for mode and current delta
// delta must be in valid range
double TrilinearTraction::Strength(int mode,double delta)
{   if(mode==1)
    {   if(delta<=uI2)
        {    if(break1is2I)
                return stress1;
            else
                return (sI2*(delta-umidI)+stress1*(uI2-delta))/(uI2-umidI);
        }
        else
            return sI2*(delIc-delta)/(delIc-uI2);
    }
    if(delta<=uII2)
    {    if(break1is2II)
            return stress2;
        else
            return (sII2*(delta-umidII)+stress2*(uII2-delta))/(uII2-umidII);
    }
    else
        return sII2*(delIIc-delta)/(delIIc-uII2);
}

// Return area under the coshesive law up to u (only used in J integral)
// Assumes ue <= u <= uc
double TrilinearTraction::WorkEnergy(int mode,double u)
{   if(mode==1)
    {   if(u<=uI2)
        {   if(break1is2I)
                 return 0.5*stress1*umidI;
            else
            {   double s2 = (sI2*(u-umidI)+stress1*(uI2-u))/(uI2-umidI);
                return 0.5*(u*stress1 + s2*(u - umidI));
            }
        }
        else
        {   double s2 = (delIc-u)*sI2/(delIc-uI2);
            return 0.5*(stress1*uI2 + sI2*(u-umidI) + s2*(u-uI2));
        }
    }
    if(u<=uII2)
    {   if(break1is2II)
            return 0.5*stress2*umidII;
        else
        {   double s2 = (sII2*(u-umidII)+stress2*(uII2-u))/(uII2-umidII);
            return 0.5*(u*stress2 + s2*(u - umidII));
        }
    }
    else
    {   double s2 = (delIIc-u)*sII2/(delIIc-uII2);
        return 0.5*(stress2*uII2 + sII2*(u-umidII) + s2*(u-uII2));
    }
}

// Return dissipated energy up to delta.
// delta must be in valid range (delta>umid)
double TrilinearTraction::DissipatedEnergy(int mode,double delta)
{   if(mode==1)
    {   if(delta>=uI2)
            return JI_1c + JI_2c*(delta-uI2)/(delIc-uI2);
        else
            return JI_1c*(delta-umidI)/(uI2-umidI);
    }
    if(delta>uII2)
        return JII_1c + JII_2c*(delta-uII2)/(delIIc-uII2);
    else
        return JII_1c*(delta-umidII)/(uII2-umidII);
}

// Get D from delta (needed for MixedModeTraction)
// delta must be in valid range
double TrilinearTraction::GetDFromDelta(int mode,double delta)
{   double sr,D = 0.;
    if(mode==1)
    {   sr = sI2/stress1;
        if(delta<=uI2)
        {   if(!break1is2I)
                D = 1. - (umidI/delta)*(sr + (1-sr)*(uI2-delta)/(uI2-umidI));
        }
        else
            D = 1. - (sr*umidI/delta)*(delIc-delta)/(delIc-uI2);
    }
    else
    {   sr = sII2/stress2;
        if(delta<=uII2)
        {   if(!break1is2II)
                D = 1. - (umidII/delta)*(sr + (1-sr)*(uII2-delta)/(uII2-umidII));
        }
        else
            D = 1. - (sr*umidII/delta)*(delIIc-delta)/(delIIc-uII2);
    }
    return D;
}

// Get delta from D (needed for MixedModeTraction)
double TrilinearTraction::GetDeltaFromD(int mode,double D)
{   double delta;
    if(mode==1)
    {   delta = umidI;
        if(D<DIbreak)
        {   if(!break1is2I)
            {   delta = umidI*(stress1*uI2-sI2*umidI)/
                            (stress1*(1.-D)*uI2 + umidI*(D*stress1-sI2));
            }
        }
        else
        {   delta = umidI*delIc*sI2/
                    (stress1*(1.-D)*(delIc-uI2) + umidI*sI2);
        }
    }
    else
    {   delta = umidII;
        if(D<DIIbreak)
        {   if(!break1is2II)
            {   delta = umidII*(stress2*uII2-sII2*umidII)/
                            (stress2*(1.-D)*uII2 + umidII*(D*stress2-sII2));
            }
        }
        else
        {   delta = umidII*delIIc*sII2/
                    (stress2*(1.-D)*(delIIc-uII2) + umidII*sII2);
        }
    }
    return delta;
}

#pragma mark TrilinearTraction::Trilinear Specific Methods

// print to output window
void TrilinearTraction::PrintTriLinearModel(const char *mode,double Jc,double sigc,double ke,double ue,
                                            double sigbreak,double ubreak,double uc) const
{
    char label[10];
    PrintProperty(PackLabel(label,"G",mode,"c"),Jc*UnitsController::Scaling(0.001),UnitsController::Label(ERR_UNITS));
    PrintProperty(PackLabel(label,"sig",mode,"1"),sigc*UnitsController::Scaling(1.e-6),UnitsController::Label(PRESSURE_UNITS));
    if(kI1>0.)
        PrintProperty(PackLabel(label,"k",mode,""),ke*UnitsController::Scaling(1.e-6),UnitsController::Label(TRACTIONSLOPE_UNITS));
    PrintProperty(PackLabel(label,"upk",mode,"1"),ue,UnitsController::Label(CULENGTH_UNITS));
    cout <<  endl;
    
    PrintProperty(PackLabel(label,"sig",mode,"2"),sigbreak*UnitsController::Scaling(1.e-6),UnitsController::Label(PRESSURE_UNITS));
    PrintProperty(PackLabel(label,"upk",mode,"2"),ubreak,UnitsController::Label(CULENGTH_UNITS));
    PrintProperty(PackLabel(label,"u",mode,"c"),uc,UnitsController::Label(CULENGTH_UNITS));
    cout <<  endl;
}

//	Calculate properties used in analyses - here triangular law
//		k1 is slope up (F/L^3)
//		u1 is displacement at s1 peak stress (relative to u3 or 0 to 1 on input, dimenensioned on output)
//		u2 is displacement at s2 stress (relative to u3 or 0 to 1 on input, dimenensioned on output)
//		u3 is final displacement (L)
//		s1 is first peak stress (F/L^2) = k1 u1 u3
//		s2 is second break point stress (F/L^2)
//		G is toughness (F/L) = (1/2) u3*(s1 u2 + s2(1 - u1))
const char *TrilinearTraction::SetTLTractionLaw(double &s1,double &k1,double &u1,double &s2,double &u2,
                                                double &u3,double &G,bool &break1is2,double &G1,double &G2)
{
	// see if initial stiffness was used and convert to possible s1 and u1
	// if k1 is provided, s1 or u1 (but not both) must be there too
    // To get u1, need u3 as well
	if(k1>0.)
	{	if(s1>=0. && u1>=0.)
			return "Can only specify two of sigmaI(II), kI(II)e, and umidI(II) for each mode";
		else if(s1<0. && u1<0.)
			return "Must specify either sigmaI(II) or umidI(II) for a mode when you specify kI(II)e";
		else if(s1<0.)
        {   // because relative, need u3 too
            if(u3<0.)
                return "Must specify delI(II)c when umidI(II) and kI(II)e are specified";
			s1=k1*u1*u3;
        }
		else
        {   if(u3<0.)
                return "Must specify delI(II)c when sigmaI(II) and kI(II)e are specified";
			u1=s1/(k1*u3);
        }
	}
    
    // if k1 not provided (<=0) it is calculated, but it cannot because one of s1 and u1 must be nonzero
	else if(s1<=0. && u1<=0.)
		return "Must specify at least one of sigmaI(II), kI(II)e, and umidI(II) for each mode";
	
	// Check not overspecified
	int nknown=0;
	if(s1>=0.) nknown++;
	if(u1>=0.) nknown++;
	if(s2>=0.) nknown++;
	if(u2>=0.) nknown++;
	if(u3>0.) nknown++;
	if(G>0.) nknown++;
	if(nknown!=5)
		return "Must specify exactly five of sigmaI(II), delpI(II), kI(II)e, sigmaI(II)2, delpkI(II)2, delI(II)c, and JI(II)c for each mode";
	
	// only one is unknown, find the other as needed, check each one
	
	if(G<0.)
		G = 0.5*u3*(s1*u2 + s2*(1.-u1));
	
	else if(s1<0.)
	{	if(u2>0.)
			s1 = (2.*G/u3 - s2*(1.-u1))/u2;
		else
			s1 = s2;		// linear softening from s2 at u=0 down to 0 at u3
	}
	
	else if(u2<0.)
	{	if(s1<=0.)
			return "When sigmaI or sigmaII is zero, you must specify delpkI2 or delpkII2, respectively";
		u2 = (2.*G/u3 - s2*(1.-u1))/s1;
	}
	
	else if(s2<0.)
	{	if(u1>=1.)
			s2=s1;			// linear to s1 than drops to zero
		else
			s2 = (2.*G/u3 - s1*u2)/(1.-u1);
	}
	
	else if(u1<0.)
	{	if(s2<=0.)
			return "When sigmaI2 or sigmaII2 is zero, you must specify delpkI or delpkII, respectively";
		u1 = -(2.*G/u3 - s1*u2 - s2)/s2;
		if(u1<0. || u1>u2)
			return "Calculated delpkI(II) is not between 0 and delpkI(II)2";
	}
	
	else if(u3<0.)
	{	if(s1*u2 + s2*(1.-u1)<=0.)
			return "Unable to calculate delIc or delIIc because input parameters define traction law with zero or negative area.";
		u3 = 2.*G/(s1*u2 + s2*(1.-u1));
	}
	
	// convert u1 and u2 to mm
	u1 *= u3;
	u2 *= u3;
	
	// check valid displacement values
	if(u3<u2 || u3<u1 || u2<u1)
	{	return "The displacement break points must be 0 ≤ delpkI ≤ delpkI2 ≤ 1 and 0 ≤ delpkII ≤ delpkII2 ≤ 1";
	}
	if(u3<0. || s1<0. || s2<0. || G<=0.)
	{	return "The critical cod or one of the two break stresses is negative or total area is less than or equal to zero";
	}
	
	// Get initial stiffness (if appropriate)
	if(u1>0.)
		k1 = s1/u1;
	else
    {   // initial linear softening no longer allowed
        k1 = -1.;
        return "Trilinear traction cannot begin with infinite initial stiffness";
    }
    
    // constant terms
    break1is2 = (bool)DbleEqual(u1,u2);
    
    // sub areas
    G2 = 0.5*s2*u3;
    G1 = G - G2;

	return NULL;
}

// Get delta'(D) from current D
// Only called by subclass MixedMode when using Newton's method
double TrilinearTraction::GetTLDeltaPrimeFromD(double D,double Dbreak,double sigma,double s2,double ue,double u2,double uc,bool break1is2)
{   double deltaprime = 0.;
    if(D<Dbreak)
    {   if(!break1is2)
        {   double denom = sigma*(1.-D)*u2 + ue*(D*sigma-s2);
            deltaprime = ue*sigma*(u2-ue)*(sigma*u2-s2*ue)/(denom*denom);
        }
    }
    else
    {   double denom = sigma*(1.-D)*(uc-u2) + ue*s2;
        deltaprime = ue*uc*s2*sigma*(uc-u2)/(denom*denom);
    }
    return deltaprime;
}

#pragma mark TrilinearTraction::Accessors

// return material type
const char *TrilinearTraction::MaterialType(void) const { return "Trilinear Cohesive Zone"; }

