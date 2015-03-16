/********************************************************************************
    CubicTraction.cpp
    nairn-mpm-fea
    
    Created by John Nairn on 4/6/08.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/CubicTraction.hpp"
#include "Cracks/CrackSegment.hpp"

extern double mtime;

#pragma mark CubicTraction::Constructors and Destructors

// Constructors with arguments 
CubicTraction::CubicTraction(char *matName) : TractionLaw(matName)
{
}

#pragma mark CubicTraction::Initialization

/* calculate properties used in analyses - here cubic traction law
	The law is sigma = (27 smax/(4 umax^3)) u (umax-u)^2 = k u (umax-u)^2
	The toughness is Jc = (9000/16) umax smax
*/
const char *CubicTraction::VerifyAndLoadProperties(int np)
{
	const char *msg=SetTractionLaw(stress1,delIc,JIc,kI1);
	if(msg!=NULL) return msg;
	
	msg=SetTractionLaw(stress2,delIIc,JIIc,kII1);
	if(msg!=NULL) return msg;
	
	return TractionLaw::VerifyAndLoadProperties(np);
}

// print to output window
void CubicTraction::PrintMechanicalProperties(void) const
{
	PrintProperty("GcI",JIc/1000.,"J/m^2");
	PrintProperty("sigI",stress1*1.e-6,"");
	PrintProperty("uI",delIc,"mm");
	PrintProperty("kI0",1.e-6*kI1*delIc*delIc,"MPa/mm");
    cout <<  endl;
	
	PrintProperty("GcII",JIIc/1000.,"J/m^2");
	PrintProperty("sigII",stress2*1.e-6,"");
	PrintProperty("uII",delIIc,"mm");
	PrintProperty("kII0",1.e-6*kII1*delIIc*delIIc,"MPa/mm");
    cout <<  endl;
	
	PrintProperty("n",nmix,"");
	cout << endl;
}

// history variables are the current peak elastic displacement in mode I or mode II
char *CubicTraction::InitHistoryData(void)
{
	double *p = CreateAndZeroDoubles(2);
    return (char *)p;
}

#pragma mark CubicTraction::Traction Law

// Traction law - assume trianglar shape with unloading from down slope back to the origin
void CubicTraction::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,double dx,double dy,double area)
{
	double Tn=0.,Tt=0.,GI=0.,GII=0.;
	double *upeak =(double *)cs->GetHistoryData();
	
	// normal force and GI (only if open)
	if(nCod>0.)
	{	// is it failed?
		if(nCod>delIc)
		{	cs->SetMatID(0);			// then debonded
			GI = JIc;
		}
		else
		{	if(nCod>upeak[0]) upeak[0]=nCod;		// new peak reached
			double keff=kI1*(delIc-upeak[0])*(delIc-upeak[0]);
			Tn=keff*nCod;
			
			// get GI for failure law
			double d=nCod/delIc;
			GI=JIc*d*d*(6.-8.*d+3.*d*d);		// N/m
		}
	}
	
	// shear force and GII always
    // is it failed?
    double absTCod=fabs(tCod);
    if(absTCod>delIIc)
    {	cs->SetMatID(0);			// then debonded
        GII = JIIc;
    }
    else
    {	if(absTCod>upeak[1]) upeak[1]=absTCod;		// new peak reached in either direction
        double keff=kII1*(delIIc-upeak[1])*(delIIc-upeak[1]);
        Tt=keff*tCod;
        
        // get GII for failure law
        double d=tCod/delIIc;
        GII=JIIc*d*d*(6.-8.*d+3.*d*d);		// N/m
    }
	
    // failure criterion
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
// return released enegery = total energy - recoverable energy (due to elastic unloading)
//		when fullEnergy is false
// units of N/mm
double CubicTraction::CrackTractionEnergy(CrackSegment *cs,double nCod,double tCod,bool fullEnergy)
{
	double tEnergy=0.;
	
	// always get entire area under the curve
	
	// normal energy only if opened
	if(nCod>0.)
	{	double d=nCod/delIc;
		tEnergy=JIc*d*d*(6.-8.*d+3.*d*d)/1000.;		// N/mm
	}
	
	// shear energy always
	double d=tCod/delIIc;
	tEnergy+=JIIc*d*d*(6.-8.*d+3.*d*d)/1000.;		// N/mm
	
	// subtract recoverable energy when want released energy
	if(!fullEnergy)
	{	double *upeak=(double *)cs->GetHistoryData();
		double keff;
		if(nCod>0.)
		{	double Tn=0.;
			keff=kI1*(delIc-upeak[0])*(delIc-upeak[0]);
			Tn=keff*nCod*1.e-6;
			tEnergy-=0.5*Tn*nCod;					// N/mm
		}
		
		// shear energy always
		double Tt=0.;
		keff=kII1*(delIIc-upeak[1])*(delIIc-upeak[1]);
		Tt=keff*tCod*1.e-6;
		tEnergy-=0.5*Tt*tCod;						// N/mm
	}

	return tEnergy;
}

#pragma mark CubicTraction::Accessors

// Return the material tag
int CubicTraction::MaterialTag(void) const { return CUBICTRACTIONMATERIAL; }

// return material type
const char *CubicTraction::MaterialType(void) const { return "Cubic Traction Law"; }

// set traction law on initialization
const char *CubicTraction::SetTractionLaw(double &smax,double &umax,double &G,double &k1)
{
	// specify smax and G
	if(umax<0.)
	{	if(smax<0. || G<0.)
			return "Too few cubic traction law properties were supplied.";
		umax=16.*G/(9000.*smax);
	}
	
	// specify umax and G
	else if(smax<0.)
	{	if(G<0.)
			return "Too few cubic traction law properties were supplied.";
		smax=16.*G/(9000.*umax);
	}
	
	// specify umax and smax
	else if(G<0.)
	{	G=9000.*umax*smax/16.;
	}
	
	// specified them all, which is an error
	else
	{	return "Must supply exactly two of delIc, sigmaI, JIc and exactly two of delIIc, sigmaII, JIIc.";
	}
	
	// stress prefactors to get force
	// Multiply by 1e6 to get N/mm^3/mm^2 (kg-m/sec^2/mm/mm^4) to g-mm/sec^2 / mm^3 / mm^2
	k1=27.*smax*1.e6/(4.*umax*umax*umax);
	
	return NULL;
}




