/********************************************************************************
    CubicTraction.cpp
    nairn-mpm-fea
 
	The Cubic law is
		
		sigma = (27 smax/(4 umax^3)) u (umax-u)^2 = k u (umax-u)^2
 
		k = (27 smax/(4 umax^3))
 
	The toughness is
 
		G = (9/16) umax smax
    
    Created by John Nairn on 4/6/08.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/CubicTraction.hpp"
#include "Cracks/CrackSegment.hpp"
#include "System/UnitsController.hpp"

extern double mtime;

#pragma mark CubicTraction::Constructors and Destructors

// Constructor
CubicTraction::CubicTraction(char *matName,int matID) : TractionLaw(matName,matID)
{
}

#pragma mark CubicTraction::Initialization

// Set each mode in the law
const char *CubicTraction::VerifyAndLoadProperties(int np)
{
	const char *msg=SetCubicTractionLaw(stress1,delIc,JIc,kI1);
	if(msg!=NULL) return msg;
	
	msg=SetCubicTractionLaw(stress2,delIIc,JIIc,kII1);
	if(msg!=NULL) return msg;
	
	return TractionLaw::VerifyAndLoadProperties(np);
}

// print to output window
void CubicTraction::PrintMechanicalProperties(void) const
{
    PrintCubicModel("I",JIc,stress1,delIc,kI1);
    PrintCubicModel("II",JIIc,stress2,delIIc,kII1);
	PrintProperty("n",nmix,"");
	cout << endl;
}

#pragma mark CubicTraction::History Data Methods

// history variables are the current peak elastic displacement in mode I or mode II
char *CubicTraction::InitHistoryData(char *pchr)
{
    numTractionHistory = 2;
	double *p = CreateAndZeroDoubles(pchr,2);
    return (char *)p;
}

#pragma mark CubicTraction::Traction Law

// Traction law - assume trianglar shape with unloading from down slope back to the origin
void CubicTraction::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,Vector *n,Vector *t,double area)
{
	double Tn=0.,Tt=0.,GI=0.,GII=0.;
	double *upeak =(double *)cs->GetHistoryData();
	
	// normal force and GI (only if open)
	if(nCod>0.)
	{	// is it failed?
		if(nCod>delIc)
		{	upeak[CUBIC_DELN]=delIc;
			cs->SetMatID(0);			// then debonded
		}
		else
		{	if(nCod>upeak[CUBIC_DELN]) upeak[CUBIC_DELN]=nCod;		// new peak reached
			double keff=kI1*(delIc-upeak[CUBIC_DELN])*(delIc-upeak[CUBIC_DELN]);
			Tn=keff*nCod;
		}
	}
	
	// shear force and GII always
    // is it failed?
    double absTCod=fabs(tCod);
    if(absTCod>delIIc)
	{	upeak[CUBIC_DELT]=delIIc;
		cs->SetMatID(0);			// then debonded
    }
    else
    {	if(absTCod>upeak[CUBIC_DELT]) upeak[CUBIC_DELT]=absTCod;		// new peak reached in either direction
        double keff=kII1*(delIIc-upeak[CUBIC_DELT])*(delIIc-upeak[CUBIC_DELT]);
        Tt=keff*tCod;
    }

    // get total dissipated energy
	CrackDissipatedEnergy(cs,GI,GII);

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
	
	// force is traction times area projected onto plane of unit vectors (units F)
	// tract = -area*(Tn*n + Tt*t)
	// In 2D, if t=(dx,dy), then n=(-dy,dx)
	cs->tract.x = -area*(Tn*n->x + Tt*t->x);
	cs->tract.y = -area*(Tn*n->y + Tt*t->y);
	cs->tract.z = -area*(Tn*n->z + Tt*t->z);
}

// Return current traction law strain energy (Int T.du).
//	This energy is needed for J integral (and only used in J Integral)
// units of F/L
double CubicTraction::CrackWorkEnergy(CrackSegment *cs,double nCod,double tCod)
{
	double d,workEnergy;
	
	// normal energy only if opened
	if(nCod>0.)
	{	d = nCod/delIc;
		workEnergy = JIc*d*d*(6. - 8.*d + 3.*d*d);
	}
	else
		workEnergy = 0.;
	
	// shear energy always
	d = fabs(tCod)/delIIc;
	workEnergy += JIIc*d*d*(6. - 8.*d + 3.*d*d);
	
	return workEnergy;
}

// Get mode I and II energy that has been released by current segment
// For this decoupled law it is given by
//     Int (1/2)phi(delta)delta
// It is used when reporting energy released during crack growth or decohesion
//   and added to Jtip to report total dissipated energy by the crack growth.
// Also used by this law to detect failure.
// units of F/L
void CubicTraction::CrackDissipatedEnergy(CrackSegment *cs,double &GI,double &GII)
{
	// Get area under curve up to maximum reached
	double *upeak = (double *)cs->GetHistoryData();
	
	// Get normal energy
	double d=upeak[CUBIC_DELN]/delIc;
	GI = JIc*(4.-3.*d)*d*d*d;
	
	// Get normal energy
	d=upeak[CUBIC_DELT]/delIIc;
	GII = JIIc*(4.-3.*d)*d*d*d;
    
    // set on segment
    cs->czmdG.x = GI;
    cs->czmdG.y = GII;
    if(upeak[CUBIC_DELN]>0 || upeak[CUBIC_DELT]>0)
        cs->czmdG.z = 1.;
}

#pragma mark CubicTraction::Accessors

// return material type
const char *CubicTraction::MaterialType(void) const { return "Cubic Traction Law"; }





