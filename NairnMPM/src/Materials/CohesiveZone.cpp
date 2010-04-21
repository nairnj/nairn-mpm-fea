/********************************************************************************
    CohesiveZone.cpp
    NairnMPM
    
    Created by John Nairn on 3/21/08.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/CohesiveZone.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include "Cracks/CrackSegment.hpp"

extern double mtime;

#pragma mark CohesiveZone::Constructors and Destructors

// Constructors with arguments 
CohesiveZone::CohesiveZone(char *matName) : TractionLaw(matName)
{
	// mode I cohesive law (all others set to -1 in superclasses)
	kI1=-1;				// kIe
	
	// mode II cohesive law (all others set to -1 in superclasses)
	kII1=-1;			// kIIe
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
	
	return SetTractionLaw(stress2,kII1,delIIc,JIIc,umidII);
}

// print to output window
void CohesiveZone::PrintMechanicalProperties(void)
{
	PrintProperty("GcI",JIc,"J/m^2");
	PrintProperty("sigI",stress1,"");
	PrintProperty("uI",delIc,"mm");
	PrintProperty("kI",1.0e-6*kI1,"MPa/mm");
	PrintProperty("upkI",umidI,"mm");
    cout <<  endl;

	PrintProperty("GcII",JIIc,"J/m^2");
	PrintProperty("sigII",stress2,"");
	PrintProperty("uII",delIIc,"mm");
	PrintProperty("kII",1.0e-6*kII1,"MPa/mm");
	PrintProperty("upkII",umidII,"mm");
    cout <<  endl;
}

// history variable is the current peak elastic displacement in mode I or mode II
// initialize to peak of triangle (meaning elastic up to the peak)
char *CohesiveZone::MaterialData(void)
{
    // allocate two doubles
    char *p=new char[2*sizeof(double)];
    double *h=(double *)p;
	
	h[0]=umidI;
	h[1]=umidII;
	
    return p;
}

#pragma mark CohesiveZone::Traction Law

// Traction law - assume trianglar shape with unloading from down slope back to the origin
void CohesiveZone::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,double dx,double dy,double area)
{
	double Tn=0.,Tt=0.;
	double *upeak =(double *)cs->GetHistoryData();
	
	// normal force (only if open)
	if(nCod>0.)
	{	// is it failed?
		if(nCod>delIc)
		{	cs->SetMatID(0);			// then debonded
			cout << "# debond at t=" << 1000.*mtime << " and (x,y) = (" << cs->x << "," <<cs-> y << ")" << endl;
		}
		else
		{	if(nCod>upeak[0]) upeak[0]=nCod;			// new peak reached
			double keff=kI1*((delIc-upeak[0])*umidI)/((delIc-umidI)*upeak[0]);
			Tn=keff*nCod;
		}
	}
	
	// shear (if still bonded)
	if(cs->MatID()>=0)
	{	// is it failed?
		if(tCod>delIIc)
		{	cs->SetMatID(0);			// then debonded
			cout << "# shear debond t=" << 1000.*mtime << " and (x,y) = (" << cs->x << "," <<cs-> y << ")" << endl;
			Tn=0.;						// turn off normal traciton when failed
		}
		else
		{	if(fabs(tCod)>upeak[1]) upeak[1]=fabs(tCod);			// new peak reached either direction
			double keff=kII1*((delIIc-upeak[1])*umidII)/((delIIc-umidII)*upeak[1]);
			Tt=keff*tCod;
		}
	}
	
	// force is traction time area projected onto x-y plane
	cs->tract.x=area*(Tn*dy - Tt*dx);
	cs->tract.y=area*(-Tn*dx - Tt*dy);
}

// return total energy (which is needed for path independent J) under traction law curve
//		when fullEnergy is true
// return released enegery = total energy - recoverable energy (due to elastic unloading)
//		when fullEnergy is false
// units of N/mm
double CohesiveZone::CrackTractionEnergy(CrackSegment *cs,double nCod,double tCod,bool fullEnergy)
{
	double tEnergy=0.;
	
	// always get entire area under the curve
	
	// normal energy only if opened
	if(nCod>0.)
	{	if(nCod<umidI)
		{	double Tn=kI1*nCod*1.e-6;				// now in units of N/mm^2
			tEnergy=0.5*Tn*nCod;					// N/mm
		}
		else
		{	double s2=(delIc-nCod)*stress1/(delIc-umidI);				// stress in N/mm^2
			tEnergy=0.5*(umidI*stress1 + (nCod-umidI)*(stress1+s2));	// N/mm
		}
	}
	
	// shear energy always
	if(fabs(tCod)<umidII)
	{	double Tt=kII1*tCod*1.e-6;					// now in units of N/mm^2
		tEnergy+=0.5*Tt*tCod;						// N/mm
	}
	else
	{	double s2=(delIIc-fabs(tCod))*stress2/(delIIc-umidI);					// stress in N/mm^2
		tEnergy+=0.5*(umidII*stress2 + (fabs(tCod)-umidII)*(stress2+s2));		// N/mm
	}
	
	// subtract recoverable energy when want released energy
	if(!fullEnergy)
	{	double *upeak=(double *)cs->GetHistoryData();
		double keff;
		if(nCod>0.)
		{	double Tn=0.;
			keff=kI1*((delIc-upeak[0])*umidI)/((delIc-umidI)*upeak[0]);
			Tn=keff*nCod*1.e-6;
			tEnergy-=0.5*Tn*nCod;					// N/mm
		}
		
		// shear energy always
		double Tt=0.;
		keff=kII1*((delIIc-upeak[1])*umidII)/((delIIc-umidII)*upeak[1]);
		Tt=keff*tCod*1.e-6;
		tEnergy-=0.5*Tt*tCod;						// N/mm
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
	
	u1 is displacement at peak stress (mm)
*/
const char *CohesiveZone::SetTractionLaw(double &s1,double &k1,double &u2,double &G,double &u1)
{
	// specify k1, s1, and G, but not u2
	if(u2<0.)
	{	if(s1<0. || G<0.)
			return "Too few cohesive law properties were supplied.";
		u2=G/(500.*s1);
	}
	
	// specify k1, u2, and G, but not s1
	else if(s1<0.)
	{	if(G<0.)
			return "Too few cohesive law properties were supplied.";
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
	
	// If k1 not not provided, set peak to u1/u2=.225926299 to best match Needelman Cubic law energy
	if(k1<0.)
	{	u1=0.225926299*u2;
		k1=s1/u1;
	}
	else
	{	// verify second slope is positive
		u1=s1/k1;
		if(u2<u1)
			return "The critical displacement is before the peak displacement. Increase JIc, JIIc, delIc, and/or delIIc.";
	}
		
		
	// Multiply by 1e6 to get N/mm/mm^2 (kg-m/sec^2/mm/mm^2) to g-mm/sec^2 / mm / mm^2
	k1*=1.0e6;
	
	return NULL;
}

