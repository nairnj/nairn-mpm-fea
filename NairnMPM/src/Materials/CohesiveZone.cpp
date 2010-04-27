/********************************************************************************
    CohesiveZone.cpp
    NairnMPM
 
	This material is based on mixed mode cohesive law with saw tooth traction
	laws in paper
 
	J. L. Hogberg, "Mixed Mode Cohesive Law," Int. J. Fract, v146, 549-559 (2006)
 
	It give total area under the cohesive law as
 
	   J = JIc sin^2 q + JIIc cos^2 q
 
    where q is mixity that ranges from 90 for mode I to 0 for mode II. It is
	defined by
  
	   tan q = (deln/delIc) / (delt/delIIc)
 
	or ratio of displacement normal and tangential to crack (deln and delt)
	delative to critical failure displacementin that mode (delIc and delIIc)
    
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

// history variables:
// h[0] is 1-w which is 1 for undamaged and slope*lampk after damage
// h[1] is lambda max
char *CohesiveZone::MaterialData(void)
{
    // allocate two doubles
    char *p=new char[2*sizeof(double)];
    double *h=(double *)p;
	
	h[0]=1.;
	h[1]=0.;
	
    return p;
}

#pragma mark CohesiveZone::Traction Law

// Traction law - assume trianglar shape with unloading from down slope back to the origin
void CohesiveZone::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,double dx,double dy,double area)
{
	// get history variables
	double *h =(double *)cs->GetHistoryData();
	double lambdaMax=h[1];
	double Sfxn=0.;
	
	// key terms
	if(nCod<0.) nCod=0.;
	double wbar=nCod/delIc;				// dimensionless opening displacement
	double wwp=nCod/umidI;				// displacement relative to opening peak
	double vbar=tCod/delIIc;			// dimensionless shear displacement
	double vvp=tCod/umidII;				// displacement relative to shear peak
	double lambda=sqrt(wbar*wbar+vbar*vbar);			// mixed mode displacement
	double lambdap=lambda/sqrt(wwp*wwp+vvp*vvp);		// mixed mode peak
	
	// new peak load
	if(lambda>lambdaMax)
	{	if(lambda>=1.0)
		{	//it is no failed
			cs->SetMatID(0);			// set it debonded
			cout << "# debond at t=" << 1000.*mtime << " and (x,y) = (" << cs->x << "," <<cs-> y << ")" << endl;
		}
		else if(lambda<lambdap)
		{	// below the peak for onset of damage (i.e., softening)
			Sfxn=1./lambdap;			// S/lambda
			h[1]=lambda;				// save lambdaMax
		}
		else
		{	// Past the peak into softening regime
			Sfxn=(1.-lambda)/(lambda*(1.-lambdap));			// S/lambda
			h[1]=lambda;									// save lambdaMax
			h[0]=Sfxn*lambdap;								// save unloading stiffness S(lmax)*lp/lmax
		}
	}
	
	// has unloaded, use previous stiffness from history variable
	else if(lambda>0.)
	{	Sfxn=h[0]/lambdap;				// S/lambda when unloaded
	}

	// Find normal and tangential forces
	double Tn=Sfxn*wbar*stress1;
	double Tt=Sfxn*vbar*stress2;
	
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
	
	// key terms
	if(nCod<0.) nCod=0.;
	double wbar=nCod/delIc;				// dimensionless opening displacement
	double wwp=nCod/umidI;					// displacement relative to opening peak
	double vbar=tCod/delIIc;				// dimensionless shear displacement
	double vvp=tCod/umidII;				// displacement relative to shear peak
	double lambda=sqrt(wbar*wbar+vbar*vbar);			// mixed mode displacement
	if(lambda<=0.) return tEnergy;
	double lambdap=lambda/sqrt(wwp*wwp+vvp*vvp);		// mixed mode peak
	
	// always normalized energy as area under S(lambda) curve divided by lambda^2
	// it is Jbar/lambda^2 (see Hogberg, Table 1)
	
	if(lambda<=lambdap)
	{	// before the peak
		tEnergy=0.5/lambdap;
	}
	else
	{	// after the peak
		tEnergy=0.5*(2. - lambda - lambdap/lambda)/(lambda*(1.-lambdap));
	}

	// subtract recoverable energy when want released energy
	if(!fullEnergy)
	{	double *h=(double *)cs->GetHistoryData();
		tEnergy-=0.5*h[0]/lambdap;
	}
	
	// convert to actual units (N/mm)
	// Eenergy = Jbar*(stress1*delIc*sin^2 q + stress2*delIIc*cos^2 q)
	tEnergy*=stress1*delIc*wbar*wbar+stress2*delIIc*vbar*vbar;
	
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
	{	u1*=u2;
		k1=s1/u1;
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

