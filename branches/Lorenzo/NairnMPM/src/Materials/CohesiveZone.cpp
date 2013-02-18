/********************************************************************************
    CohesiveZone.cpp
    NairnMPM
 
	This material is based on mixed mode cohesive law with saw tooth traction
	laws. For mixed mode behavior, it uses method advocated by Thouless, which
    is simple method, but appears more effective than methods that try to be more
    coupled.
 
	Basically, mode I and mode II are decoupled for calculation of GI and GII. If
    either exceeds its critical COD, that it debonds. If neither does, the criterion
 
    (GI/GIc)^n + (GII/GIIc)^n = 1
 
    is used to decide if it fails.
    
    Created by John Nairn on 3/21/08.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/CohesiveZone.hpp"
#include "Cracks/CrackSegment.hpp"

extern double mtime;

#pragma mark CohesiveZone::Constructors and Destructors

// Constructors with arguments 
CohesiveZone::CohesiveZone(char *matName) : TractionLaw(matName)
{
	// mode I cohesive law (all others set to -1 in superclasses)
	// others are: stress1,delIc,JIc
	kI1=-1.;			// kIe
	umidI=-1.;
	
	// mode II cohesive law (all others set to -1 in superclasses)
	// others are: stress2,delIIc,JIIc
	kII1=-1;			// kIIe
	umidII=-1.;
}

#pragma mark CohesiveZone::Initialization

// no properties to read
char *CohesiveZone::InputMat(char *xName,int &input)
{   
    if(strcmp(xName,"kIe")==0)
	{	input=DOUBLE_NUM;
		return((char *)&kI1);
	}
	
	else if(strcmp(xName,"kIIe")==0)
	{	input=DOUBLE_NUM;
		return((char *)&kII1);
	}
	
	else if(strcmp(xName,"delpkI")==0)
	{	input=DOUBLE_NUM;
		return((char *)&umidI);
	}
	
	else if(strcmp(xName,"delpkII")==0)
	{	input=DOUBLE_NUM;
		return((char *)&umidII);
	}
	
    return TractionLaw::InputMat(xName,input);
}

/* calculate properties used in analyses - here triangular law
	In terms of J (J/m^2) and stress (MPa)
	    umax = J/(500*stress), k = 1000 stress^2/J
	In terms of k and umax
	    J = 250 k umax^2,   stress = k umax/2
*/
const char *CohesiveZone::VerifyProperties(int np)
{
	const char *msg=SetTractionLaw(stress1,kI1,delIc,JIc,umidI);
	if(msg!=NULL) return msg;
	
	msg=SetTractionLaw(stress2,kII1,delIIc,JIIc,umidII);
	if(msg!=NULL) return msg;
	
	return TractionLaw::VerifyProperties(np);
}

// print to output window
void CohesiveZone::PrintMechanicalProperties(void)
{
	PrintProperty("GcI",JIc,"J/m^2");
	PrintProperty("sigI",stress1,"");
	PrintProperty("uIc",delIc,"mm");
	if(kI1>0.) PrintProperty("kI",1.0e-6*kI1,"MPa/mm");
	PrintProperty("upkI",umidI,"mm");
    cout <<  endl;

	PrintProperty("GcII",JIIc,"J/m^2");
	PrintProperty("sigII",stress2,"");
	PrintProperty("uIIc",delIIc,"mm");
	if(kII1>0.) PrintProperty("kII",1.0e-6*kII1,"MPa/mm");
	PrintProperty("upkII",umidII,"mm");
    cout <<  endl;
	
	PrintProperty("n",nmix,"");
	cout << endl;
	
	// Multiply by 1e6 to get N/mm/mm^2 (kg-m/sec^2/mm/mm^2) to g-mm/sec^2 / mm / mm^2
	sIc=stress1*1.e6;
	sIIc=stress2*1.e6;
}

// history variables:
// h[0] is max mode I opening (starting at first peak location)
// h[1] is max mode II opening (starting at first peak location)
char *CohesiveZone::InitHistoryData(void)
{
    double *p = CreateAndZeroDoubles(2);
	p[0]=umidI;
	p[1]=umidII;
    return (char *)p;
}

#pragma mark CohesiveZone::Traction Law

// Traction law - assume trianglar shape with unloading from down slope back to the origin
void CohesiveZone::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,double dx,double dy,double area)
{
	double Tn=0.,Tt=0.,GI=0.,GII=0.;
	double *upeak =(double *)cs->GetHistoryData();
	
	// normal force (only if open)
	if(nCod>0.)
	{	// is it failed?
		if(nCod>delIc)
		{	cs->SetMatID(0);                        // then debonded
			ReportDebond(mtime,cs,1.0);
		}
		else
		{	if(nCod>upeak[0]) upeak[0]=nCod;                        // new peak reached
			double keff=sIc*(delIc-upeak[0])/((delIc-umidI)*upeak[0]);
			Tn=keff*nCod;
			
			// get GI for failure law
			if(nCod<umidI)
				GI=0.0005*kI1*nCod*nCod;                             // now in units of N/m
			else
			{	double s2=(delIc-nCod)*stress1/(delIc-umidI);
				GI=500.*(umidI*stress1 + (nCod-umidI)*(stress1+s2));        // now in units of N/m
			}
		}
	}
	
	// shear (if still bonded)
	if(cs->MatID()>=0)
	{	// is it failed?
		double absTCod=fabs(tCod);
		if(absTCod>delIIc)
		{	cs->SetMatID(0);                        // then debonded
			ReportDebond(mtime,cs,0.0);
			Tn=0.;                                          // turn off normal traction when failed
		}
		else if(absTCod>0.)
		{	if(absTCod>upeak[1]) upeak[1]=absTCod;                    // new peak reached either direction
			double keff=sIIc*(delIIc-upeak[1])/((delIIc-umidII)*upeak[1]);
			Tt=keff*tCod;
			
			// shear energy always
			if(absTCod<umidII)
				GII=0.0005*kII1*tCod*tCod;               // now in units of N/m
			else
			{	double s2=(delIIc-absTCod)*stress2/(delIIc-umidII);
				GII=500.*(umidII*stress2 + (absTCod-umidII)*(stress2+s2));      // now in units of N/m
			}
		}
	}
	
	// mixed mode failure? (nmix<=0 uses infinity which means fails when either COD becomes critical)
	if(cs->MatID()>=0 && nmix>0)
	{	if(pow(GI/JIc,nmix)+pow(GII/JIIc,nmix) > 1)
		{	cs->SetMatID(0);				// now debonded
			ReportDebond(mtime,cs,GI/(GI+GII));
			Tn=0.;
			Tt=0.;
		}
	}
	
	// force is traction times area projected onto x-y plane
	cs->tract.x=area*(Tn*dy - Tt*dx);
	cs->tract.y=area*(-Tn*dx - Tt*dy);
}

// return total energy (which is needed for path independent J) under traction law curve
//		when fullEnergy is true
// return released energy = total energy - recoverable energy (due to elastic unloading)
//		when fullEnergy is false
// units of N/mm
double CohesiveZone::CrackTractionEnergy(CrackSegment *cs,double nCod,double tCod,bool fullEnergy)
{
	double tEnergy=0.;
	
	// always get entire area under the curve
	
	// normal energy only if opened
	if(nCod>0.)
	{	if(nCod<umidI)
		{	double Tn=kI1*nCod;
			tEnergy=0.5e-6*Tn*nCod;					// now in units of N/mm
		}
		else
		{	double s2=(delIc-nCod)*stress1/(delIc-umidI);                   // stress in N/mm^2
			tEnergy=0.5*(umidI*stress1 + (nCod-umidI)*(stress1+s2));		// now in units of N/mm
		}
	}
	
	// shear energy always
	if(fabs(tCod)<umidII)
	{	double Tt=kII1*tCod;
		tEnergy+=0.5e-6*Tt*tCod;                         // now in units of N/mm
	}
	else
	{	double s2=(delIIc-fabs(tCod))*stress2/(delIIc-umidII);                           // stress in N/mm^2
		tEnergy+=0.5*(umidII*stress2 + (fabs(tCod)-umidII)*(stress2+s2));               // now in units of N/mm
	}
	
	// subtract recoverable energy when want released energy
	if(!fullEnergy)
	{	double *upeak=(double *)cs->GetHistoryData();
		double keff;
		if(nCod>0.)
		{	keff=sIc*(delIc-upeak[0])/((delIc-umidI)*upeak[0]);
			double Tn=keff*nCod;
			tEnergy-=0.5e-6*Tn*nCod;                                // now in units of N/mm
		}
		
		// shear energy always
		keff=sIIc*(delIIc-upeak[1])/((delIIc-umidII)*upeak[1]);
		double Tt=keff*tCod;
		tEnergy-=0.5e-6*Tt*tCod;									// N/mm
	}
	
	return tEnergy;
}

#pragma mark CohesiveZone::Accessors

// return material type
const char *CohesiveZone::MaterialType(void) { return "Triangular Cohesive Zone"; }

// Return the material tag
int CohesiveZone::MaterialTag(void) { return COHESIVEZONEMATERIAL; }

/* calculate properties used in analyses - here triangular law
	k1 is slope up (MPa/mm)
	u2 is final displacement (mm)
	s1 is peak stress (MPa) = k1 u1
	G is toughness (J/m^2) = 500 s1 u2 = 500 k1 u1 u2
	u1 is displacement at peak stress (as fraction of u2)
*/
const char *CohesiveZone::SetTractionLaw(double &s1,double &k1,double &u2,double &G,double &u1)
{
	// specify k1, s1, and G, but not u2
	if(u2<0.)
	{	if(s1<0. || G<0.)
			return "Must supply exactly two of delIc, sigmaI, JIc and exactly two of delIIc, sigmaII, JIIc.";
		u2=G/(500.*s1);
	}
	
	// specify k1, u2, and G, but not s1
	else if(s1<0.)
	{	if(G<0.)
			return "Must supply exactly two of delIc, sigmaI, JIc and exactly two of delIIc, sigmaII, JIIc.";
		s1=G/(500.*u2);
	}
	
	// specify k1, u2, and s1, but not G
	else if(G<0.)
	{	G=500.*s1*u2;
	}
	
	// specified them all, which is an error
	else
	{	return "Must supply exactly two of delIc, sigmaI, JIc and exactly two of delIIc, sigmaII, JIIc.";
	}
	
	// If neither k1 nor u1 provided, set peak to u1/u2=.225926299 to best match Needelman Cubic law energy
	if(k1<0. && u1<0.)
	{	u1=0.225926299*u2;
		k1=s1/u1;
	}
	
	// if only u1 is provided (and it is assumed to be relative to u1)
	else if(k1<0.)
	{	if(u1>0.)
		{	u1*=u2;
			k1=s1/u1;
		}
		else
		{	// linear softening
			u1=0.;
			k1=-1.;
		}
	}
	
	// if k1 is provided or (both are provided)
	else if(u1<0.)
	{	u1=s1/k1;
	}
	
	// specified both, which is an error
	else
	{	return "Must supply at most one of kIe and delpkI and one of kIIe and delpkII";
	}
	
	// verify final u1 is acceptable
	if(u2<u1)
		return "The critical displacement is less than peak displacement. Increase JIc, JIIc, delIc, and/or delIIc.";
		
	// Multiply by 1e6 to get N/mm/mm^2 (kg-m/sec^2/mm/mm^2) to g-mm/sec^2 / mm / mm^2
	k1*=1.0e6;
	
	return NULL;
}

