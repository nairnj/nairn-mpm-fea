/********************************************************************************
	CoulombFriction.cpp
	nairn-mpm-fea

	Friction slicing or stick contact

	Created by John Nairn, Oct 24, 3015.
	Copyright (c) 2015 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Materials/CoulombFriction.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"

#pragma mark CoulombFriction::Constructors and Destructors

// Constructors
CoulombFriction::CoulombFriction() {}

CoulombFriction::CoulombFriction(char *matName) : ContactLaw(matName)
{
	frictionCoeff = 0.0;			// <0 is stick
	frictionCoeffStatic = -1.;		// ignored if negative or if > frictionCoeff or if frictionCoeff < 0
}

#pragma mark CoulombFriction::Initialization

// Read material properties by name (in xName). Set input to variable type
// (DOUBLE_NUM or INT_NUM) and return pointer to the class variable
// (cast as a char *)
char *CoulombFriction::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    if(strcmp(xName,"coeff")==0)
	{	input=DOUBLE_NUM;
		return (char *)&frictionCoeff;
	}
	
    if(strcmp(xName,"coeffStatic")==0)
	{	input=DOUBLE_NUM;
		return (char *)&frictionCoeffStatic;
	}
	
	// does not all any from MaterialBase
    return ContactLaw::InputMaterialProperty(xName,input,gScaling);
}

// Verify input properties do calculations; if problem return string with an error message
// Don't pass on to material base
const char *CoulombFriction::VerifyAndLoadProperties(int np)
{
	if(frictionCoeff<0.)
	{	frictionStyle = STICK;
		frictionCoeffStatic = -1.;
	}
	
	else if(DbleEqual(frictionCoeff,(double)0.))
		frictionStyle = FRICTIONLESS;
	
	else
		frictionStyle = FRICTIONAL;
	
	// if>0, must be less than dynamic coefficient and convert to frictional is dynamic coeff is zero
	if(frictionCoeffStatic>0.)
	{	if(frictionCoeffStatic<frictionCoeff)
			return "The static coefficient of friction cannot be lower than the sliding coefficient of friction";
		if(frictionStyle == FRICTIONLESS)
			frictionStyle = FRICTIONAL;
	}
	
	// must call super class
	return ContactLaw::VerifyAndLoadProperties(np);
}

// print contact law details to output window
void CoulombFriction::PrintContactLaw(void) const
{
	char hline[200];
	
	switch(frictionStyle)
	{	case STICK:
			cout << "Contact by stick conditions, separation is free" << endl;
			break;
		
		case FRICTIONLESS:
			cout << "Contact by frictionless sliding" << endl;;
			break;
		
		default:
			sprintf(hline,"Contact by Coulomb friction with coefficient of friction: %.6f",frictionCoeff);
			cout << hline << endl;
			break;
	}
	
	if(frictionCoeffStatic>0.)
	{	sprintf(hline,"                      and static coefficient of friction: %.6f",frictionCoeffStatic);
		cout << hline << endl;
	}
}

#pragma mark CoulombFriction:Step Methods

// Adjust change in momentum for frictional contact in tangential direction
// If has component of tangential motion, calculate force depending on whether it is sticking or sliding
// When frictional sliding, find tangential force (times dt) and set flag, if not set flag false
// When friction heating is on, set Ftdt term and set hasFriction to true
//     (hasFriction (meaning has frictional heating value) must be initialized to false when called)
// contactArea only provided if frictional law needs it
// Normally in contact when called, but some laws might want call even when not in contact. If return
//		value is false, the contact should be treated as no contact
bool CoulombFriction::GetFrictionalDeltaMomentum(Vector *delPi,Vector *norm,double dotn,double *mredDE,double mred,
												 bool getHeating,double contactArea,bool inContact,double deltime,Vector *at) const
{
	// indicate no frictional heat yet
	*mredDE=-1.;
	
	// stick and frictionless are easy and no heating
	if(frictionStyle==STICK)
	{	// stick conditions no change
		return true;
	}
	
	else if(frictionStyle==FRICTIONLESS)
	{	// remove tangential term
		CopyScaleVector(delPi,norm,dotn);
		return true;
	}
	
	// Rest implements friction sliding
	// The initial delPi = (-N A dt) norm + (Sstick A dt) tang = dotn norm + dott tang
	// where N is normal traction (positive in compression), A is contact area, and dt is timestep
		
    // get unnormalized tangential vector and its magnitude
	// tang = delPi - dotn norm
    Vector tang;
    CopyVector(&tang,delPi);
    AddScaledVector(&tang,norm,-dotn);
    double tangMag = sqrt(DotVectors(&tang,&tang));
    
    // if has tangential motion, we need to change momemtum if frictional sliding is occuring
    if(!DbleEqual(tangMag,0.))
    {	ScaleVector(&tang,1./tangMag);
        double dott = DotVectors(delPi,&tang);
        
        // make it positive for comparison to the positive frictional force Sslide
        if(dott < 0.)
        {	ScaleVector(&tang,-1.);
            dott = -dott;
        }
		
		// Let frictional sliding force be Sslide Ac dt = f(N) Ac dt
		// Then if dott > Sslide Ac dt (which means Sstick>Sslide)
		//	a. Set delPi = dotn norm + Sslide Ac dt tang
		//      For example, Coulomb friction has Fslide dt = mu(dotn) so delPi = (norm - mu tang) dotn
		double SslideAcDt = GetSslideAcDt(-dotn,dott,0.,mred,contactArea,inContact,deltime);
		if(!inContact) return false;
		
		if(dott > SslideAcDt)
		{	CopyScaleVector(delPi,norm,dotn);
			AddScaledVector(delPi,&tang,SslideAcDt);
			
            // get frictional heating term as friction work times reduced mass
            // As heat source need Energy/sec or divide by timestep*reduced mass
            // Note: only add frictional heating during momentum update (when friction
            //   force is appropriate) and only if frictional contact heat is enabled.
            if(getHeating)
			{	if(at!=NULL)
				{	double Vs = SslideAcDt*(dott-SslideAcDt);
					AddScaledVector(at, delPi, -1./deltime);
					double AsDt = SslideAcDt*DotVectors(at,&tang)*deltime;
					if(AsDt>Vs)
					{	//*mredDE = Vs*(Vs/AsDt-1.)+0.5*AsDt;
						*mredDE = 0.5*Vs*Vs/AsDt;
					}
					else
						*mredDE = Vs - 0.5*AsDt;
				}
				else
					*mredDE = SslideAcDt*(dott-SslideAcDt);
           }
		}
    }
	
	// still in contact
	return true;
}

// Return Sslide Ac dt = f(N) Ac dt
// Input is N Ac dt (and is always positive when in contact)
// If needed in the friction law, Ac is the contact area (not yet provided)
// The relative sliding speed after correcting the momentum will be (SStickAcDt-SslideAcDt)/mred
double CoulombFriction::GetSslideAcDt(double NAcDt,double SStickAcDt,double Ac,double mred,
									  double contactArea,bool &inContact,double deltime) const
{
	// check static coefficient
	if(frictionCoeffStatic>0.)
	{	double Sstatic = frictionCoeffStatic*NAcDt;
		if(SStickAcDt<Sstatic) return Sstatic;			// it will stick
	}
	
	// S = mu N so S Ac dt = mu N Ac dt
	return frictionCoeff*NAcDt;
}

#pragma mark CoulombFriction::Accessors

// return unique, short name for this material
const char *CoulombFriction::MaterialType(void) const { return "Coulomb Friction"; }

// Return the material tag
int CoulombFriction::MaterialTag(void) const { return COULOMBFRICTIONLAW; }

// Set coefficient of friction
void CoulombFriction::SetFrictionCoeff(double newCoeff) { frictionCoeff = newCoeff; }

// True to ignore contact and revert to single velocity field
bool CoulombFriction::IgnoreContact(void) const { return false; }


