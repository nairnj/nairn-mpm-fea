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

#include "Materials/TrilinearTraction.hpp"
#include "Cracks/CrackSegment.hpp"

extern double mtime;

#pragma mark TrilinearTraction::Constructors and Destructors

// Constructors with arguments 
TrilinearTraction::TrilinearTraction(char *matName) : CohesiveZone(matName)
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
char *TrilinearTraction::InputMat(char *xName,int &input)
{
	if(strcmp(xName,"sigmaI2")==0)
	{	input=DOUBLE_NUM;
		return((char *)&sI2);
	}
	
	else if(strcmp(xName,"sigmaII2")==0)
	{	input=DOUBLE_NUM;
		return((char *)&sII2);
	}
	
	else if(strcmp(xName,"delpkI2")==0)
	{	input=DOUBLE_NUM;
		return((char *)&uI2);
	}
	
	else if(strcmp(xName,"delpkII2")==0)
	{	input=DOUBLE_NUM;
		return((char *)&uII2);
	}
	
    return CohesiveZone::InputMat(xName,input);
}

// Calculate properties used in analyses - here trilinear law
const char *TrilinearTraction::VerifyAndLoadProperties(int np)
{
	const char *msg=SetTLTractionLaw(stress1,kI1,umidI,sI2,uI2,delIc,JIc);
	if(msg!=NULL) return msg;
	
	msg=SetTLTractionLaw(stress2,kII1,umidII,sII2,uII2,delIIc,JIIc);
	if(msg!=NULL) return msg;
	
	const char *err = TractionLaw::VerifyAndLoadProperties(np);
	
	// Multiply by 1e6 to get N/mm/mm^2 (kg-m/sec^2/mm/mm^2) to g-mm/sec^2 / mm / mm^2
	sIc=stress1*1.e6;
	sIIc=stress2*1.e6;
	sIc2=sI2*1.e6;
	sIIc2=sII2*1.e6;
	break1is2I = DbleEqual(umidI,uI2);
	break1is2II = DbleEqual(umidII,uII2);
	
	return err;
}

// print to output window
void TrilinearTraction::PrintMechanicalProperties(void) const
{
	PrintProperty("GcI",JIc,"J/m^2");
	PrintProperty("sigI1",stress1,"");
	if(kI1>0.) PrintProperty("kI",1.0e-6*kI1,"MPa/mm");
	PrintProperty("upkI1",umidI,"mm");
    cout <<  endl;
	PrintProperty("sigI2",sI2,"");
	PrintProperty("upkI2",uI2,"mm");
	PrintProperty("uIc",delIc,"mm");
    cout <<  endl;
	
	PrintProperty("GcII",JIIc,"J/m^2");
	PrintProperty("sigII1",stress2,"");
	if(kII1>0.) PrintProperty("kII",1.0e-6*kII1,"MPa/mm");
	PrintProperty("upkII1",umidII,"mm");
    cout <<  endl;
	PrintProperty("sigII2",sII2,"");
	PrintProperty("upkII2",uII2,"mm");
	PrintProperty("uIIc",delIIc,"mm");
    cout <<  endl;
	
	PrintProperty("n",nmix,"");
	cout << endl;
}

#pragma mark TrilinearTraction::Traction Law

// Traction law - assume trianglar shape with unloading from down slope back to the origin
void TrilinearTraction::CrackTractionLaw(CrackSegment *cs,double nCod,double tCod,double dx,double dy,double area)
{
	double Tn=0.,Tt=0.,GI=0.,GII=0.;
	double *upeak =(double *)cs->GetHistoryData();
	
	// normal force and GI (only if open)
	if(nCod>0.)
	{	// is it failed?
		if(nCod>delIc)
		{	cs->SetMatID(0);                        // then debonded
			GI = JIc;
		}
		else
		{	if(nCod>upeak[0]) upeak[0]=nCod;                        // new peak reached
			double keff;
			if(upeak[0]<=uI2)
			{	if(break1is2I)
					keff = sIc/upeak[0];
				else
					keff = (sIc2*(upeak[0]-umidI)+sIc*(uI2-upeak[0]))/((uI2-umidI)*upeak[0]);
			}
			else
				keff=sIc2*(delIc-upeak[0])/((delIc-uI2)*upeak[0]);
			Tn=keff*nCod;
			
			// get GI for failure law
			if(nCod<umidI)
            {   // note that initial linear softening never because umidI=0
				GI=0.0005*kI1*nCod*nCod;                             // now in units of N/m
            }
			else if(nCod<=uI2)
			{	if(break1is2I)
					GI=500.*stress1*umidI;
				else
				{	double s2=(sI2*(nCod-umidI)+stress1*(uI2-nCod))/(uI2-umidI);
					GI=500.*(nCod*stress1 + s2*(nCod - umidI));
				}
			}
			else
			{	double s2=(delIc-nCod)*sI2/(delIc-uI2);
				GI=500.*(stress1*uI2 + sI2*(nCod-umidI) + s2*(nCod-uI2));        // now in units of N/m
			}
		}
	}
	
	// shear force and GII always
    // is it failed?
    double absTCod=fabs(tCod);
    if(absTCod>delIIc)
    {	cs->SetMatID(0);                        // then debonded
        GII = JIIc;
    }
    else if(absTCod>0.)
    {	if(absTCod>upeak[1]) upeak[1]=absTCod;                        // new peak reached
        double keff;
        if(upeak[1]<=uII2)
        {	if(break1is2II)
                keff = sIIc/upeak[1];
            else
                keff = (sIIc2*(upeak[1]-umidII)+sIIc*(uII2-upeak[1]))/((uII2-umidII)*upeak[1]);
        }
        else
            keff=sIIc2*(delIIc-upeak[1])/((delIIc-uII2)*upeak[1]);
        Tt=keff*tCod;
        
        // get GII for failure law
        if(absTCod<umidII)
        {   // note that initial linear softening never because umidI=0
            GII=0.0005*kII1*tCod*tCod;                             // now in units of N/m
        }
        else if(absTCod<=uII2)
        {	if(break1is2II)
                GII=500.*stress2*umidII;
            else
            {	double s2=(sII2*(absTCod-umidII)+stress2*(uII2-absTCod))/(uII2-umidII);
                GII=500.*(absTCod*stress2 + s2*(absTCod - umidII));
            }
        }
        else
        {	double s2=(delIIc-absTCod)*sII2/(delIIc-uII2);
            GII=500.*(stress2*uII2 + sII2*(absTCod-umidII) + s2*(absTCod-uII2));        // now in units of N/m
        }
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
			Tn=0.;                                       // turn off in tractions, if calculated
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
double TrilinearTraction::CrackTractionEnergy(CrackSegment *cs,double nCod,double tCod,bool fullEnergy)
{
	double tEnergy=0.;
	
	// always get entire area under the curve
	
	// normal energy only if opened
	if(nCod>0.)
	{	if(nCod<umidI)
        {   // note that initial linear softening never because umidI=0
			double Tn=kI1*nCod;
			tEnergy=0.5e-6*Tn*nCod;					// now in units of N/mm
		}
		else if(nCod<=uI2)
		{	if(break1is2II)
				tEnergy=0.5*stress1*umidII;
			else
			{	double s2=(sI2*(nCod-umidI)+stress1*(uI2-nCod))/(uI2-umidI);                // stress in N/mm^2
				tEnergy=0.5*(nCod*stress1 + s2*(nCod - umidI));									// now in units of N/mm
			}
		}
		else
		{	double s2=(delIc-nCod)*sI2/(delIc-uI2);
			tEnergy=0.5*(stress1*uI2 + sI2*(nCod-umidI) + s2*(nCod-uI2));        // now in units of N/m,
		}
	}
	
	// shear energy always
	double absTCod=fabs(tCod);
	if(absTCod<umidII)
    {   // note that initial linear softening never because umidI=0
		double Tt=kII1*tCod;
		tEnergy+=0.5e-6*Tt*tCod;                         // now in units of N/mm
	}
	else if(absTCod<=uII2)
	{	if(break1is2II)
			tEnergy+=0.5*stress2*umidII;
		else
		{	double s2=(sII2*(absTCod-umidII)+stress2*(uII2-absTCod))/(uII2-umidII);
			tEnergy+=0.5*(absTCod*stress2 + s2*(absTCod - umidII));						 // now in units of N/mm
		}
	}
	else
	{	double s2=(delIIc-absTCod)*sII2/(delIIc-uII2);
		tEnergy+=0.5*(stress2*uII2 + sII2*(absTCod-umidII) + s2*(absTCod-uII2));        // now in units of N/mm
	}
	
	// subtract recoverable energy when want released energy
	if(!fullEnergy)
	{	double *upeak=(double *)cs->GetHistoryData();
		double keff;
		if(nCod>0.)
		{	if(upeak[0]<=uI2)
			{	if(break1is2I)
					keff = sIc/upeak[0];
				else
					keff = (sIc2*(upeak[0]-umidI)+sIc*(uI2-upeak[0]))/((uI2-umidI)*upeak[0]);
			}
			else
				keff=sIc2*(delIc-upeak[0])/((delIc-uI2)*upeak[0]);
			double Tn=keff*nCod;
			tEnergy-=0.5e-6*Tn*nCod;                                // now in units of N/mm
		}
		
		// shear energy always
		if(upeak[1]<=uII2)
		{	if(break1is2II)
				keff = sIIc/upeak[1];
			else
				keff = (sIIc2*(upeak[1]-umidII)+sIIc*(uII2-upeak[1]))/((uII2-umidII)*upeak[1]);
		}
		else
			keff=sIIc2*(delIIc-upeak[1])/((delIIc-uII2)*upeak[1]);
		double Tt=keff*tCod;
		tEnergy-=0.5e-6*Tt*tCod;									// N/mm
	}
	
	return tEnergy;
}

#pragma mark TrilinearTraction::Accessors

// return material type
const char *TrilinearTraction::MaterialType(void) const { return "Trilinear Cohesive Zone"; }

// Return the material tag
int TrilinearTraction::MaterialTag(void) const { return TRILINEARTRACTIONMATERIAL; }

/* Calculate properties used in analyses - here triangular law
	k1 is slope up s1 (MPa/mm)
	s1 is peak stress (MPa) = k1 u1
	u1 is displacement at s1 peak stress (relative to u3 or 0 to 1)
	s2 is second break point stress (MPa)
	u2 is displacement at s2 stress (relative to u3 or 0 to 1)
	u3 is final displacement (mm)
	G is toughness (J/m^2) = 500 (s1 u2 _ s2 u3 - s2 u1)
*/
const char *TrilinearTraction::SetTLTractionLaw(double &s1,double &k1,double &u1,double &s2,double &u2,double &u3,double &G)
{
	// see if initial stiffness was used and convert to possible s1 and u1
	// if k1 is provided, s1 or u1 (but not both) must be there too
	if(k1>0.)
	{	if(s1>=0. && u1>=0.)
			return "Can only specify two of sigmaI(II), kI(II)e, and umidI(II) for each mode";
		else if(s1<0. && u1<0.)
			return "Must specify either sigmaI(II) or umidI(II) for a mode where you specify kI(II)e";
		else if(s1<0.)
			s1=k1*u1;
		else
			u1=s1/k1;
	}
    
    // if k1 not provided (<1) it is calculated, but it cannot but one of s1 and u1 must be nonzero
    // to start with linear softening, need u1=0 and s1>0
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
		G = 500.*u3*(s1*u2 + s2*(1.-u1));
	
	else if(s1<0.)
	{	if(u2>0.)
			s1 = (0.002*G/u3 - s2*(1.-u1))/u2;
		else
			s1 = s2;		// linear softening from s2 at u=0 down to 0 at u3
	}
	
	else if(u2<0.)
	{	if(s1<=0.)
			return "When sigmaI or sigmaII is zero, you must specify delpkI2 or delpkII2, respectively";
		u2 = (0.002*G/u3 - s2*(1.-u1))/s1;
	}
	
	else if(s2<0.)
	{	if(u1>=1.)
			s2=s1;			// linear to s1 than drops to zero
		else
			s2 = (0.002*G/u3 - s1*u2)/(1.-u1);
	}
	
	else if(u1<0.)
	{	if(s2<=0.)
			return "When sigmaI2 or sigmaII2 is zero, you must specify delpkI or delpkII, respectively";
		u1 = -(0.002*G/u3 - s1*u2 - s2)/s2;
		if(u1<0. || u1>u2)
			return "Calculated delpkI(II) is not between 0 and delpkI(II)2";
	}
	
	else if(u3<0.)
	{	if(s1*u2 + s2*(1.-u1)<=0.)
			return "Unable to calculate delIc or delIIc because input parameters define traction law with zero or negative area.";
		u3 = 0.002*G/(s1*u2 + s2*(1.-u1));
	}
	
	// convert u1 and u2 to mm
	u1*=u3;
	u2*=u3;
	
	// check valid displacement values
	if(u3<u2 || u3<u1 || u2<u1)
	{	return "The displacement break points must be 0 ≤ delpkI ≤ delpkI2 ≤ 1 and 0 ≤ delpkII ≤ delpkII2 ≤ 1";
	}
	if(u3<0. || s1<0. || s2<0. || G<=0.)
	{	return "The critical cod or one of the two break stresses is negative or total area is less than or equal to zero";
	}
	
	// Get initial stiffness (if appropriate)
	// Multiply by 1e6 to get N/mm/mm^2 (kg-m/sec^2/mm/mm^2) to g-mm/sec^2 / mm / mm^2
	if(u1>0.)
		k1 = s1/u1;
	else
		k1=-1.;
	k1*=1.0e6;
	
	return NULL;
}

