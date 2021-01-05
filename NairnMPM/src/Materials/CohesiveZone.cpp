/********************************************************************************
    CohesiveZone.cpp
    nairn-mpm-fea
 
    This material is parent class of any traction law that has an initial
    linear elastic region. It adds parameters for slope and displacement
    at end of linear elatic regime. This instands implements a sawtooth
    (or triangular law).
 
    This material uses decoupled mode I and mode II (e.g., see work of
    Thouless and others). If either mode exceeds its critical COD, that it debonds.
    If neither does, this decoupled approach need a failure criterion. The
    one used is
 
    (GI/GIc)^n + (GII/GIIc)^n = 1
 
    Seom subclasses of this class implemeted coupled methods
 
    Created by John Nairn on 3/21/08.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/CohesiveZone.hpp"
#include "Materials/CoupledSawTooth.hpp"
#include "Cracks/CrackSegment.hpp"
#include "System/UnitsController.hpp"

extern double mtime;

#pragma mark CohesiveZone::Constructors and Destructors

// Constructor
CohesiveZone::CohesiveZone(char *matName,int matID) : TractionLaw(matName,matID)
{
	// mode I cohesive law (all others set to -1 in superclasses)
	// others are: stress1,delIc,JIc
	kI1=-1.;			// initial elastic slope (keep <0 for linear softening)
	umidI=-1.;			// peak
	
	// mode II cohesive law (all others set to -1 in superclasses)
	// others are: stress2,delIIc,JIIc
	kII1=-1;			// initial elastic slope (keep <0 for linear softening)
	umidII=-1.;			// peak
}

#pragma mark CohesiveZone::Initialization

// Read properteis (read read in super classes)
char *CohesiveZone::InputTractionLawProperty(char *xName,int &input,double &gScaling)
{   
    if(strcmp(xName,"kIe")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&kI1,gScaling,1.e6);
	}
	
	else if(strcmp(xName,"kIIe")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&kII1,gScaling,1.e6);
	}
	
	else if(strcmp(xName,"delpkI")==0)
	{	input=DOUBLE_NUM;
		return((char *)&umidI);
	}
	
	else if(strcmp(xName,"delpkII")==0)
	{	input=DOUBLE_NUM;
		return((char *)&umidII);
	}
	
    return TractionLaw::InputTractionLawProperty(xName,input,gScaling);
}

// Calculate properties used in analyses - here triangular law
// Do mode I and mode II separately
const char *CohesiveZone::VerifyAndLoadProperties(int np)
{
	const char *msg=SetSawToothTractionLaw(stress1,kI1,delIc,JIc,umidI,phiI1,RI1);
	if(msg!=NULL) return msg;
	
	msg=SetSawToothTractionLaw(stress2,kII1,delIIc,JIIc,umidII,phiII1,RII1);
	if(msg!=NULL) return msg;
	
	// go to parent
	return TractionLaw::VerifyAndLoadProperties(np);
}

// print to output window
void CohesiveZone::PrintMechanicalProperties(void) const
{
    PrintSawToothModel("I",JIc,stress1,delIc,kI1,umidI,-1.);
    PrintSawToothModel("II",JIIc,stress2,delIIc,kII1,umidII,-1.);
	if(MaterialID()!=COUPLEDSAWTOOTHMATERIAL)
    {   PrintProperty("n",nmix,"");
        cout << endl;
    }
}

#pragma mark CohesiveZone::History Data Methods

// history variables:
// h[0=CZ_DELN] is max mode I opening (starting at first peak location)
// h[1=CZ_DELT] is max mode II opening (starting at first peak location)
char *CohesiveZone::InitHistoryData(char *pchr)
{
    numTractionHistory = 2;
    double *p = CreateAndZeroDoubles(pchr,2);
	p[CZ_DELN]=umidI;
	p[CZ_DELT]=umidII;
    return (char *)p;
}

#pragma mark CohesiveZone::Basic Functions

// Return the strength for mode and current delta
// delta must be in valid range
double CohesiveZone::Strength(int mode,double delta)
{   if(mode==1)
        return stress1*(delIc-delta)/(delIc-umidI);
    return stress2*(delIIc-delta)/(delIIc-umidII);
}

// Return area under the cohesive law up to u (only used in J integral)
// Assumes ue <= u <= uc
double CohesiveZone::WorkEnergy(int mode,double u)
{   if(mode==1)
        return 0.5*stress1*(u*(delIc-u) + delIc*(u-umidI))/(delIc-umidI);
    return 0.5*stress2*(u*(delIIc-u) + delIIc*(u-umidII))/(delIIc-umidII);
}

// Return dissipated energy up to delta.
// delta must be in valid range (delta>umid)
double CohesiveZone::DissipatedEnergy(int mode,double delta)
{   if(mode==1)
        return JIc*(delta-umidI)/(delIc-umidI);
    return JIIc*(delta-umidII)/(delIIc-umidII);
}

// Return the derivative of strength with repect to delta
// delta must be in valid range
double CohesiveZone::StrengthPrime(int mode,double delta)
{   if(mode==1)
        return -stress1/(delIc-umidI);
    return -stress2/(delIIc-umidII);
}

// Get D from delta (needed for MixedModeTraction)
// delta must be in valid range
double CohesiveZone::GetDFromDelta(int mode,double delta)
{   double delbar,up;
    if(mode==1)
    {   delbar = delta/delIc;
        up = umidI/delIc;
    }
    else
    {   delbar = delta/delIIc;
        up = umidII/delIIc;
    }
    return 1. - (up/(1.-up))*((1-delbar)/delbar);
}

// Get delta from D (needed for MixedModeTraction)
double CohesiveZone::GetDeltaFromD(int mode,double D)
{   double up;
    if(mode==1)
    {   up = umidI/delIc;
        return delIc*up/(1.-D*(1.-up));
    }
    up = umidII/delIIc;
    return delIIc*up/(1.-D*(1.-up));
}

#pragma mark CohesiveZone::Traction Law

// Traction law - This method should be general for any cohesive law that has an
// initial elastic regime provided delta parameters are initialized to end of the elastic regime
// This class implements sawtooth law
void CohesiveZone::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,Vector *n,Vector *t,double area)
{
	double Tn=0.,Tt=0.,GI=0.,GII=0.;
	double *upeak = (double *)cs->GetHistoryData();
	
	// normal force and GI (only if open)
	if(nCod>0.)
	{	// is it failed?
		if(nCod>delIc)
		{	upeak[CZ_DELN]=delIc;
			cs->SetMatID(0);                        // then debonded
		}
		else
		{	if(nCod>upeak[CZ_DELN]) upeak[CZ_DELN]=nCod;                        // new peak reached
            double keff = Strength(1,upeak[CZ_DELN])/upeak[CZ_DELN];
			Tn = keff*nCod;
        }
	}
	
	// shear force and GII always
    // is it failed?
    double absTCod = fabs(tCod);
    if(absTCod>delIIc)
	{	upeak[CZ_DELT] = delIIc;
		cs->SetMatID(0);                                // then debonded
    }
    else if(absTCod>0.)
    {	if(absTCod>upeak[CZ_DELT]) upeak[CZ_DELT] = absTCod;          // new peak reached either direction
        double keff = Strength(2,upeak[CZ_DELT])/upeak[CZ_DELT];
        Tt=keff*tCod;
    }
	
    // get total dissipated energy
	CrackDissipatedEnergy(cs,GI,GII);
	
    if(cs->MatID()<0)
    {   // it failed above in pure mode
        ReportDebond(mtime, cs, GI/(GI+GII),GI+GII);
        Tn = 0.;                                       // turn off in tractions, if calculated
        Tt = 0.;
    }
	else if(nmix>0)
    {   // mixed mode failure? (nmix<=0 uses infinity which means fails when either COD becomes critical)
		if(pow(GI/JIc,nmix)+pow(GII/JIIc,nmix) > 1)
		{	cs->SetMatID(0);				// now debonded
			ReportDebond(mtime,cs,GI/(GI+GII),GI+GII);
			Tn = 0.;                                       // turn off in tractions, if calculated
			Tt = 0.;
		}
	}
	
	// force is traction times area projected onto plane of unit vectors (units F)
	// tract = -area*(Tn*n + Tt*t)
	// In 2D, if t=(dx,dy), then n=(-dy,dx)
	cs->tract.x = -area*(Tn*n->x + Tt*t->x);
	cs->tract.y = -area*(Tn*n->y + Tt*t->y);
	cs->tract.z = -area*(Tn*n->z + Tt*t->z);
}

// Return current traction law work energy (Int T.du) up to current u=(nCod,tCod)
// Only called by J unless subclass failure to override DissipatedEnergy()
// This method should be general for any cohesive law that has an initial elastic regime
// provided delta parameters are initialized to end of the elastic regime.
// This class implements sawtooth law
double CohesiveZone::CrackWorkEnergy(CrackSegment *cs,double nCod,double tCod)
{
	double workEnergy;
	
	// normal energy only if opened
	if(nCod>0.)
	{	if(nCod<umidI || delIc==umidI)
		{   // elastic regime
			double Tn = kI1*nCod;
			workEnergy = 0.5*Tn*nCod;
		}
		else
        {   workEnergy = WorkEnergy(1,nCod);
        }
	}
	else
		workEnergy = 0.;
	
	// add shear energy always
	double absTCod = fabs(tCod);
	if(absTCod<umidII || delIIc==umidII)
	{   // elastic loading
		double Tt = kII1*tCod;
		workEnergy += 0.5*Tt*tCod;
	}
	else
    {   workEnergy += WorkEnergy(2,absTCod);
    }
	
	return workEnergy;
}

// Get mode I and II energy that has been released by current segment
// For this decoupled law it is given by Int (1/2)phi(delta)delta
// It is used when reporting energy released during crack growth or decohesion
//      and added to Jtip to report total dissipated energy by the crack growth.
// Also used by this law to detect failure.
// This should be general for any law that has an initial elastic regime provided
//      deltas are initialized to end of the elastic regime.
// This class implements an sawtooth law
void CohesiveZone::CrackDissipatedEnergy(CrackSegment *cs,double &GI,double &GII)
{
	// Get area under curve up to maximum reached
	double *upeak = (double *)cs->GetHistoryData();
	
    // Get normal energy (only if damage)
    if(upeak[CZ_DELN]>umidI)
	{	GI = DissipatedEnergy(1,upeak[CZ_DELN]);
        cs->czmdG.z = 1.;       // flag for disspated energy
	}
    else
        GI = 0.;
	
	// Get shear energy (only if damage)
	if(upeak[CZ_DELT]>umidII)
	{	GII = DissipatedEnergy(2,upeak[CZ_DELT]);
        cs->czmdG.z = 1.;       // flag for dissipated energy
	}
    else
        GII = 0.;
    
    // add to crack properties
    cs->czmdG.x = GI;
    cs->czmdG.y = GII;
}

#pragma mark CohesiveZone::SawTooth Specific (no overrides)

//	Calculate properties used in analyses - here triangular law
//		k1 is slope up (F/L^3)
//		u2 is final displacement (L)
//		u1 is displacement at peak stress (as fraction of u2 on input, but dimensioned on output)
//		s1 is peak stress (F/L^2) = k1 u1 (output u1)
//		G is toughness (F/L) = (1/2) s1 u2 = (1/2) k1 u1 u2 (output u1)
//      phi1 is constant terms for varphi function
//      R1 is constant terms for R(delta) function
const char *CohesiveZone::SetSawToothTractionLaw(double &s1,double &k1,double &u2,double &G,
                                                 double &u1,double &phi1,double &R1)
{
	// specify s1 and G, but not u2
	if(u2<0.)
	{	if(s1<0. || G<0.)
			return "Must supply exactly two of delIc, sigmaI, JIc and exactly two of delIIc, sigmaII, JIIc.";
		u2 = 2.*G/s1;
	}
	
	// specify u2 and G, but not s1
	else if(s1<0.)
	{	if(G<0.)
			return "Must supply exactly two of delIc, sigmaI, JIc and exactly two of delIIc, sigmaII, JIIc.";
		s1 = 2.*G/u2;
	}
	
	// specify u2 and s1, but not G
	else if(G<0.)
	{	G = 0.5*s1*u2;
	}
	
	// specified them all, which is an error
	else
	{	return "Can only supply exactly two of delIc, sigmaI, JIc and exactly two of delIIc, sigmaII, JIIc.";
	}
	
	// If neither k1 nor u1 provided, set peak to u1/u2=.225926299 to best match Needelman Cubic law energy
	if(k1<0. && u1<0.)
	{	u1 = 0.225926299*u2;
		k1 = s1/u1;
	}
	
	// if only u1 is provided (and it is assumed to be relative to u2 on input)
	else if(k1<0.)
	{	if(u1>0.)
		{	u1 *= u2;
			k1 = s1/u1;
		}
		else
		{	// linear softening (no longer allowed
			u1 = 0.;
			k1 = -1.;
            return "Triangular traction material does not allow infinite initial stiffness";
		}
	}
	
	// if k1 is provided alone
	else if(u1<0.)
    {   if(k1>0.)
            u1 = s1/k1;
        else
            return "Triangular traction material does not allow zero initial stiffness";
	}
	
	// specified both, which is an error
	else
	{	return "Can only supply one of kIe and delpkI and one of kIIe and delpkII";
	}
	
	// verify final u1 is acceptable
	if(u2<u1)
		return "The critical displacement is less than peak displacement. Increase JIc, JIIc, delIc, and/or delIIc.";
    
    // constant terms
    phi1 = s1*u2/(u2-u1);
    R1 = u1/(u2-u1);

	return NULL;
}

// print to output window
// alpha only applies to exponential law
void CohesiveZone::PrintSawToothModel(const char *mode,double Jc,double sigc,
                                      double uc,double ke,double ue,double alpha) const
{
    char label[10];
    PrintProperty(PackLabel(label,"G",mode,"c"),Jc*UnitsController::Scaling(0.001),UnitsController::Label(ERR_UNITS));
    PrintProperty(PackLabel(label,"sig",mode,""),sigc*UnitsController::Scaling(1.e-6),
                  UnitsController::Label(PRESSURE_UNITS));
    PrintProperty(PackLabel(label,"u",mode,"c"),uc,UnitsController::Label(CULENGTH_UNITS));
    cout <<  endl;
    // note that ke<=0 or infinite stiffness is not longer allowed
    if(ke>0.)
    {   PrintProperty(PackLabel(label,"k",mode,""),ke*UnitsController::Scaling(1.e-6),
                      UnitsController::Label(TRACTIONSLOPE_UNITS));
    }
    PrintProperty(PackLabel(label,"upk",mode,""),ue,UnitsController::Label(CULENGTH_UNITS));
    if(alpha>0.)
        PrintProperty(PackLabel(label,"alpha",mode,""),alpha,"");
    cout <<  endl;
}

// Get delta' for SawTooth law at D
// Only called by subclass MixedMode traction law if using Newton's method
double CohesiveZone::GetSawToothDeltaPrimeFromD(double D,double uc,double ue)
{   double up = ue/uc;
    double denom = 1.-D*(1.-up);
    return uc*up*(1.-up)/(denom*denom);
}

#pragma mark CohesiveZone::Exponential Specifics (here so available to Coupled and Mixed Mode)

//    Calculate properties used in analyses - here exponential law
//        G is toughness
//        u2 is final displacement (L)
//              exactly one of above two required
//        k1 is slope up (F/L^3)
//        u1 is displacement at peak stress (as fraction of u2 on input, but dimensioned on output)
//        s1 is peak stress (F/L^2) = k1 u1 (output u1)
//              exactly two of above three required
//        fealpha = 2(1/alpha - exp(-alpha)/(1-exp(-alpha)) is exponential factor
//              required between 0 and 1
const char *CohesiveZone::SetExponentialTractionLaw(double &s1,double &k1,double &u2,double &G,double &u1,double &fealpha)
{
    // Must have G or u2, but not both
    if((G>0. && u2>0.) || (G<=0. && u2<=0.))
        return "Exponential cohesive laws must set Jc or delc, and never both";
    
    // must have two of s1, k1, u1
    int style=0;
    if(s1>0. && k1>0. && u1<=0.)
        style = 1;      // has s1 and k1
    else if(s1>0. && k1<=0. && u1>0.)
        style = 2;      // has s1 and u1
    else if(s1<=0. && k1>0. && u1>0.)
        style = 3;      // has k1 and u1
    if(style==0)
        return "Exponential cohesive laws must set exactly two of ke, sigma, and delpk";

    // provide uc
    if(u2>0.)
    {   switch(style)
        {   case 1:
                // provided s1 and k1
                u1 = s1/(k1*u2);
                break;
            case 2:
                // provided s1 and u1
                k1 = s1/(u1*u2);
                break;
            default:
                // provided k1 and u1
                s1 = k1*u1*u2;
                break;
        }
        
        // need toughness
        G = 0.5*s1*u2*(u1 + (1.-u1)*fealpha);
    }
    
    // provided Jc
    else
    {   switch(style)
        {   case 1:
                // provided s1 and k1
                u2 = (2.*k1*G - (1.-fealpha)*s1*s1)/(k1*s1*fealpha);
                u1 = s1/(k1*u2);
                break;
            case 2:
                // provided s1 and u1
                u2 = 2.*G/(s1*(u1 + (1.-u1)*fealpha));
                k1 = s1/(u1*u2);
                break;
            default:
                // provided k1 and u1
                u2 = sqrt(2.*G/(k1*u1*(u1 + (1.-u1)*fealpha)));
                s1 = k1*u1*u2;
                break;
        }
    }
    
    // convert u1 to dimesions
    u1 *= u2;

    return NULL;
}


#pragma mark CohesiveZone::Accessors

// return material type
const char *CohesiveZone::MaterialType(void) const { return "Triangular Cohesive Zone"; }

