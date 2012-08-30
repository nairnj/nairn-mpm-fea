/********************************************************************************
    CoupledSawTooth.hpp
    NairnMPM

    Created by John Nairn on 8/30/12.
    Copyright (c) 2102 John A. Nairn, All rights reserved.
 
    This material is based on mixed-mode damage law such that normal and shear
    tractions are
 
       Tn = (1-D) k dn      and Tt = (1-D) k dt
 
    where D = (dc (deff - dpk)) / (deff (dc - dpk)) is a damge parameter a
    D = 0 if deff < dpk. Here deff = |(dn,dt)|, dc is critical deff at
    failure, and dpk is deff at the peak of Teff vs deff where
    Teff = |(Tm,Tt)| and Teff = (1-D) k deff
 
    If we require (Tn, Tt) to be gradient of a potential (which this law does)
    we must have single k above (same stiffness in normal and traction). This law
    implies pure mode I and mode II tractions are idential (same cohesive stress
    of k dpk and same critical COD of dc).
 
    Failure occurs when G = (1/2) k dpk dc, independent of mode and cohesive stress
    if sc = k dpk
 
    Created by John Nairn on 8/30/12.
    Copyright (c) 2102 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/CoupledSawTooth.hpp"
#include "Cracks/CrackSegment.hpp"

extern double mtime;

#pragma mark CoupledSawTooth::Constructors and Destructors

// Constructors with arguments 
CoupledSawTooth::CoupledSawTooth(char *matName) : CohesiveZone(matName)
{
	// use mode I settings of superclass law
}

#pragma mark CoupledSawTooth::Initialization

/* calculate properties used in analyses - here triangular law
    In terms of J (J/m^2) and stress (MPa)
    umax = J/(500*stress), k = 1000 stress^2/J
    In terms of k and umax
    J = 250 k umax^2,   stress = k umax/2
*/
const char *CoupledSawTooth::VerifyProperties(int np)
{
	const char *msg=SetTractionLaw(stress1,kI1,delIc,JIc,umidI);
	if(msg!=NULL) return msg;
	
    // do not need to call base material class methods
	return NULL;
}

// print to output window
void CoupledSawTooth::PrintMechanicalProperties(void)
{
	PrintProperty("Gc",JIc,"J/m^2");
	PrintProperty("sig",stress1,"");
	PrintProperty("uc",delIc,"mm");
	if(kI1>0.) PrintProperty("k",1.0e-6*kI1,"MPa/mm");
	PrintProperty("upk",umidI,"mm");
    cout <<  endl;
    
	// Multiply by 1e6 to get N/mm/mm^2 (kg-m/sec^2/mm/mm^2) to g-mm/sec^2 / mm / mm^2
	sIc=stress1*1.e6;
}

// history variables:
// h is max effective displacement opening (starting at peak location)
char *CoupledSawTooth::MaterialData(void)
{
    double *h=new double;
    *h=0.;
    return (char *)h;
}

#pragma mark CohesiveZone::Traction Law

// Traction law - assume trianglar shape with unloading from down slope back to the origin
void CoupledSawTooth::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,double dx,double dy,double area)
{
	double Tn=0.,Tt=0.;
	double *upeak =(double *)cs->GetHistoryData();
    double deff = sqrt(nCod*nCod + tCod*tCod);
    
    // is it debonded
    if(deff > delIc)
    {   cs->SetMatID(0);                        // now debonded
        ReportDebond(mtime,cs,1.0);
        cs->tract.x = 0.;
        cs->tract.y = 0.;
        return;
    }
    
    // is it a new peak?
    if(deff > umidI) upeak[0] = deff;
    
    // stiffness same for both modes keff = (1-D)k = (1-D) sc (df-d)/(d*(df-d0)
    // Note: prior to deff reaching d0, upeak[0]=umidI and  keff = sc/d0 = k
    double keff=sIc*(delIc-upeak[0])/((delIc-umidI)*upeak[0]);
    
    // normal force (only if open, closed handled by crack contact)
    if(nCod>0.) Tn = keff*nCod;
    
    Tt = keff*tCod;
	
	// force is traction times area projected onto x-y plane
	cs->tract.x=area*(Tn*dy - Tt*dx);
	cs->tract.y=area*(-Tn*dx - Tt*dy);
}

// return total energy (which is needed for path independent J) under traction law curve
//		when fullEnergy is true
// return released energy = total energy - recoverable energy (due to elastic unloading)
//		when fullEnergy is false
// units of N/mm
double CoupledSawTooth::CrackTractionEnergy(CrackSegment *cs,double nCod,double tCod,bool fullEnergy)
{
	double tEnergy=0.;
    
    double deff = sqrt(nCod*nCod + tCod*tCod);
    
    if(deff < umidI)
    {
        
    }
    
    double G = sIc*(deff*deff - 2.*deff*delIc + umidI*delIc)/(2.*(umidI - delIc));
		
	return tEnergy;
}

#pragma mark CohesiveZone::Accessors

// return material type
const char *CoupledSawTooth::MaterialType(void) { return "Coupled Triangular Law"; }

// Return the material tag
int CoupledSawTooth::MaterialTag(void) { return COUPLEDSAWTOOTHMATERIAL; }

