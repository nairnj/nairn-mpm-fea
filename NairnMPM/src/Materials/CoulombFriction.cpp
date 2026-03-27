/********************************************************************************
	CoulombFriction.cpp
	nairn-mpm-fea

	Friction sliding or stick contact

	Created by John Nairn, Oct 24, 3015.
	Copyright (c) 2015 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/CoulombFriction.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "System/UnitsController.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"

#pragma mark CoulombFriction::Constructors and Destructors

// Constructor
CoulombFriction::CoulombFriction(char *matName,int matID) : ContactLaw(matName,matID)
{
	frictionCoeff = 0.0;			// <0 is stick
	frictionCoeffStatic = -1.;		// ignored if negative or if frictionCoeff < 0
	
	// lasrt used in Revision 3491. If entered printing message it is ignored
	displacementOnly = 0.;
	Dc = -1.;
}

#pragma mark CoulombFriction::Initialization

// Read material properties by name (in xName). Set input to variable type
// (DOUBLE_NUM or INT_NUM) and return pointer to the class variable
// (cast as a char *)
char *CoulombFriction::InputContactProperty(char *xName,int &input,double &gScaling)
{
    if(strcmp(xName,"coeff")==0)
	{	input=DOUBLE_NUM;
		return (char *)&frictionCoeff;
	}
	
    else if(strcmp(xName,"coeffStatic")==0)
	{	input=DOUBLE_NUM;
		return (char *)&frictionCoeffStatic;
	}
	
	else if(strcmp(xName,"displacementOnly")==0)
	{	// ignored - warning if enter other than 0
		input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&displacementOnly,gScaling,1.e6);
	}
	
	else if(strcmp(xName,"Dc")==0)
	{	// ignored - warning if positive
		input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&Dc,gScaling,1.e6);
	}
	
	// does not all any from MaterialBase
    return ContactLaw::InputContactProperty(xName,input,gScaling);
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
	size_t hsize=200;
	
	switch(frictionStyle)
	{	case STICK:
			cout << "Contact by stick conditions, separation is free" << endl;
			break;
		
		case FRICTIONLESS:
			cout << "Contact by frictionless sliding" << endl;;
			break;
		
		default:
			snprintf(hline,hsize,"Contact by Coulomb friction with coefficient of friction: %.6f",frictionCoeff);
			cout << hline << endl;
			break;
	}
	
	if(frictionCoeffStatic>0.)
	{	snprintf(hline,hsize,"                      and static coefficient of friction: %.6f",frictionCoeffStatic);
		cout << hline << endl;
	}
	
	cout << "   Detection by negative separation and stress < 0" << endl;
	
	// features last used in revision 3941
	if(displacementOnly>0.1 || displacementOnly<0.)
		cout << "    (displacementOnly option entered - no longer used)" << endl;
	if(Dc>=0.)
		cout << "    (Dc option entered - no longer used)" << endl;
}

#pragma mark CoulombFriction:Step Methods

// Change input momentum for stick (in delPi) to reflect frictional sliding contact law
//		(this called both my material contact and crack surface contact
// Input parameters (besides delPi)
//		norm = normal vector from material i to j
//		dotn = delPi.norm (precalculated)
//		deltaDotn = initial normal opening cod precalculated (=delta.norm)
//		mred = reduced mass
//		getHeating = true to calculate frictional heating term
//		contactArea = contact area, which is only needed by some laws
//		deltime = time step
//		delFi = force changed needed in post update calculations (only non-NULL in UPDATE_MOMENTUM_CALL)
//				and needed when frictional heating is activated
//		forCracks = true when called for crack contact
// Output
//		delPi change to reflect contact law
//		true is returned or false if decide now not in contact
//		*mredDelWf set to heat energy (actually mred*heat energy) (only if getHeating is true)
bool CoulombFriction::GetFrictionalDeltaMomentum(Vector *delPi,Vector *norm,double dotn,double deltaDotn,
							double *mredDelWf,double mred,bool getHeating,double contactArea,
							double deltime,Vector *delFi,NodalPoint *ndptr,bool forCracks) const
{
	// indicate no frictional heat yet
	*mredDelWf=-1.;
	
	// stick and frictionless are easy and no heating
	if(frictionStyle==STICK)
	{	// stick conditions no change to delPi
		return true;
	}

	// Adjust deltaDotn after mommentum update (small change, maybe not needed and
	//    it is skipped for linear or logistic regressions beause they are based
	//    on two closest particles and not total nodal force, but always done for cracks)
	if(delFi!=NULL && (forCracks || !mpmgrid.UsingRegressionForContact()))
	{	// dotn/mred = delta v, fn/mred = delta a so we want
		// extra disp = delta v*dt - 1/2 delta a dt^2 or in current terms as here:
		double fn = DotVectors(delFi,norm);
		deltaDotn += deltime*(dotn - 0.5*fn*deltime)/mred;
	}

	// provisional setting for in contact
	// Both negative separation (deltaDotn<0) and compression (dotn<0)
	bool inContact = deltaDotn<0. && dotn<0. ? true : false ;
	
	// Handle frictionless is special case. It uses provional contact
	// and adjust delPi when in contact. Frictionles has no heating.
	if(frictionStyle==FRICTIONLESS)
	{	// contact requires negative separation and positive pressure
		if(inContact)
		{	// frictionless contact - return normal component
			CopyScaleVector(delPi,norm,dotn);
			return true;
		}

		// separated so no contact
		return false;
	}

	// Rest implements for frictional sliding
	// The initial delPi = (-N A dt) norm + (Sstick A dt) tang = dotn norm + dott tang
	// where N is normal traction (positive in compression), A is contact area, and dt is timestep
	
	// we are done if this law never changes momentum when not in contact (adhesion can cause non-free separation)
	if(!inContact && HasFreeSeparation()) return false;
	
	// Get force to stick  in tangential motion
	double dott = 0.;

    // Tangengial stick from tangential direction unit vector
	// tang||tang|| = delPi - dotn norm
    Vector tang;
    CopyVector(&tang,delPi);
    AddScaledVector(&tang,norm,-dotn);
    double tangMag = sqrt(DotVectors(&tang,&tang));
	if(tangMag>0.)
    {	ScaleVector(&tang,1./tangMag);
        dott = DotVectors(delPi,&tang);
        if(dott < 0.)
		{	// make it positive for comparison to the positive frictional force Sslide
        	ScaleVector(&tang,-1.);
            dott = -dott;
        }
	}
	
	// Get frictional sliding force to be Sslide Ac dt = f(N) Ac dt where NAcDt = -fnaDt
	double SslideAcDt = GetSslideAcDt(-dotn,dott,mred,contactArea,inContact,deltime);
	
	// If provisional inContact was false and adhesive contaact did not change it
	// then exit
	if(!inContact) return false;
	
	// if dott > Sslide Ac dt (which means Sstick>Sslide), then sliding, otherwise stick
	// but Sslide Ac dt <=0 means effectively frictionless
	if(SslideAcDt<=0.)
	{	// Normal stick condition
		CopyScaleVector(delPi,norm,dotn);
		//ZeroVector(delPi);		// fixes FricionL2 when displacement only because reverting to stress&displacement
	}
	else if(dott > SslideAcDt)
	{	// Normal stick condition
		CopyScaleVector(delPi,norm,dotn);

		// frictional terms
		// revision 3941 and older, only added in momentum update (delFi!=NULL)
		AddScaledVector(delPi,&tang,SslideAcDt);
		
		// get frictional heating term as friction work times reduced mass
		// As heat source need Del Wf/sec or divided by timestep*reduced mass
		// Note: only add frictional heating during momentum update (when friction
		//   force is appropriate) and only if frictional contact heat is enabled.
		if(getHeating)
		{	// delFi must be provided when getHeating is true (acceleration is delFi/mred, which is applied after return)
			// Vs alone is first order heating (i.e., *mredDelWf = Vs is first order method)
			double Vs = SslideAcDt*(dott-SslideAcDt);
			if(delFi!=NULL)
			{	AddScaledVector(delFi, delPi, -1./deltime);
				double AsDt = SslideAcDt*DotVectors(delFi,&tang)*deltime;
				if(AsDt>Vs)
					*mredDelWf = 0.5*Vs*Vs/AsDt;
				else
					*mredDelWf = Vs - 0.5*AsDt;
			}
			else
				*mredDelWf = Vs;
			//*mredDelWf = Vs;							// revert to first order heating
		}
	}
	else
	{	// frictional stick - leave delPi as is and no neating
	}
	
	// still in contact
	return true;
}

// Return Sslide Ac dt = f(N) Ac dt
// Input is N Ac dt (and is always positive when in contact)
double CoulombFriction::GetSslideAcDt(double NAcDt,double SStickAcDt,double mred,
									  double contactArea,bool &inContact,double deltime) const
{
	// check static coefficient
	if(frictionCoeffStatic>0.)
	{	// If force to stick is less that static coefficient, then it will stick
		double Sstatic = frictionCoeffStatic*NAcDt;
		
		// returns force >= stick so contact law will pick static force
		if(SStickAcDt<=Sstatic) return Sstatic;
		
		// it did not stick, so fall through to sliding
	}
	
	// S = mu N   or   S Ac dt = mu N Ac dt
	return frictionCoeff*NAcDt;
}

#pragma mark CoulombFriction::Accessors

// return unique, short name for this material
const char *CoulombFriction::MaterialType(void) const { return "Coulomb Friction"; }

// Set coefficient of friction
void CoulombFriction::SetFrictionCoeff(double newCoeff) { frictionCoeff = newCoeff; }

// This and subclasses handle contact
bool CoulombFriction::IgnoreContact(void) const { return false; }

// Return true is frictionless contact and no adhesion
bool CoulombFriction::IsFrictionless(void) const { return frictionStyle==FRICTIONLESS; }

// Return true is stick contact
bool CoulombFriction::IsStick(void) const { return frictionStyle==STICK; }

// If no momentum change whenever separated, then return true
// When true, GetSlideAcDt() always has inContact=true
bool CoulombFriction::HasFreeSeparation(void) const { return true; }

