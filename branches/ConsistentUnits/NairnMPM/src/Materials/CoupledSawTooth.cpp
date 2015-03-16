/********************************************************************************
    CoupledSawTooth.hpp
    nairn-mpm-fea

    Created by John Nairn on 8/30/12.
    Copyright (c) 2102 John A. Nairn, All rights reserved.
 
    This material is based on mixed-mode damage law such that normal and shear
    tractions are
 
       Tn = (1-D) k dn      and Tt = (1-D) k dt
 
    where D = (dc (deff - dpk)) / (deff (dc - dpk)) is a damage parameter and
    D = 0 if deff < dpk. Here deff = |(dn,dt)|, dc is critical deff at
    failure, and dpk is deff at the peak of Teff vs deff where
 
        Teff = |(Tm,Tt)| and Teff = (1-D) k deff
 
    In damage variable, deff is actually the maximum deff attained in the
    loading history.
 
    If we require (Tn, Tt) to be gradient of a potential (which this law does)
    we must have single k above (same stiffness in normal and traction). This law
    implies pure mode I and mode II tractions are idential (same cohesive stress
    of k dpk and same critical COD of dc).
 
    But I think potential arguments need to be clarified when damage is
    occurring.
 
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
const char *CoupledSawTooth::VerifyAndLoadProperties(int np)
{
    // set off mode I settings
	const char *msg=SetTractionLaw(stress1,kI1,delIc,JIc,umidI);
	if(msg!=NULL) return msg;
	
    // mode II not allowed
    if(stress2>0. || kII1>0. || delIIc>0. || JIIc>0. || umidII>0.)
        return "Mode II properties not allowed in Coupled Triangular Traction law.";
    
    // do not need to call base material class methods
	return NULL;
}

// print to output window
void CoupledSawTooth::PrintMechanicalProperties(void) const
{
	PrintProperty("Gc",JIc/1000.,"J/m^2");
	PrintProperty("sig",stress1*1.e-6,"");
	PrintProperty("uc",delIc,"mm");
	if(kI1>0.) PrintProperty("k",1.0e-6*kI1,"MPa/mm");
	PrintProperty("upk",umidI,"mm");
    cout <<  endl;
}

// history variables:
// h is max effective displacement opening (starting at peak location)
char *CoupledSawTooth::InitHistoryData(void)
{
	double *p = CreateAndZeroDoubles(5);
    p[0] = umidI;
    return (char *)p;
}

#pragma mark CohesiveZone::Traction Law

// Traction law - assume trianglar shape with unloading down slope back to the origin
void CoupledSawTooth::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,double dx,double dy,double area)
{
	double Tn=0.,Tt=0.;
	double *upeak =(double *)cs->GetHistoryData();
    double deff = sqrt(nCod*nCod + tCod*tCod);
    
    // is it debonded
    if(deff > delIc)
    {   cs->SetMatID(0);                        // now debonded
		
		// calculate mode mixity
        ReportDebond(mtime,cs,upeak[3]/(upeak[3]+upeak[4]),1.e-3*(upeak[3]+upeak[4]));
        cs->tract.x = 0.;
        cs->tract.y = 0.;
        return;
    }
    
    // is it a new peak? (note that upeak = max(max opening, umidI or peak stress opening)
    if(deff > upeak[0]) upeak[0] = deff;
    
    // skip if deff=zero since tractions are zero, and would cause problem if pure linear softening law when deff=0
    // (deff>0 implies upeak[0]>0 even when umidI=0 for pure linear softening)
    if(deff > 0.)
    {   // stiffness same for both modes keff = (1-D)k = sc(df-dmax)/(dmax*(df-d0)) = k d0*(df-dmax)/(dmax*(df-d0))
        // Note: prior to deff reaching d0, all dmax=upeak[0]=umidI=d0 are equal and keff = sc/d0 = k
        double keff=stress1*(delIc-upeak[0])/((delIc-umidI)*upeak[0]);
        
        // normal force (only if open, closed handled by crack contact)
        if(nCod>0.) Tn = keff*nCod;
        
        Tt = keff*tCod;
    }
	
	// track mode mixity
	// Units tracked GI (in [3]) and GII (in [4]) are microN/mm
	// Multiply by 1e-3 to get J/m^2
	if(nCod>upeak[1])
	{	upeak[3] += Tn*(nCod-upeak[1]);
		upeak[1] = nCod;
	}
	double abstCod=fabs(tCod);
	if(abstCod>upeak[2])
	{	double ddt = abstCod - upeak[2];		// d(delta_t)
		upeak[4] += fabs(Tt)*ddt;
		upeak[2] = abstCod;
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
// For this law it is from area under T vs deff curve
double CoupledSawTooth::CrackTractionEnergy(CrackSegment *cs,double nCod,double tCod,bool fullEnergy)
{
	double tEnergy=0.;
    
    double deff = sqrt(nCod*nCod + tCod*tCod);
    
    if(deff < umidI)
    {   // note that linear softening never here because umidI=0
    	double T=kI1*deff;
        tEnergy=0.5e-6*T*deff;					// now in units of N/mm
    }
    else
    {   // G = stress1*(deff*deff - 2.*deff*delIc + umidI*delIc)/(2.*(umidI - delIc));
    	double s2=(delIc-deff)*stress1/(delIc-umidI);                   // stress in N/mm^2
        tEnergy=0.5*(umidI*stress1 + (deff-umidI)*(stress1+s2));		// now in units of N/mm
    }
   
	// subtract recoverable energy when want released energy
	if(!fullEnergy && deff>0.)
	{	double *upeak =(double *)cs->GetHistoryData();
		double keff = stress1*(delIc-upeak[0])/((delIc-umidI)*upeak[0]);
        double T=keff*deff;
        tEnergy-=0.5e-6*T*deff;                                // now in units of N/mm
    }
    
	return tEnergy;
}

#pragma mark CohesiveZone::Accessors

// return material type
const char *CoupledSawTooth::MaterialType(void) const { return "Coupled Triangular Law"; }

// Return the material tag
int CoupledSawTooth::MaterialTag(void) const { return COUPLEDSAWTOOTHMATERIAL; }

