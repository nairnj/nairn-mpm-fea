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

#pragma mark CoulombFriction::Constructors and Destructors

// Constructor
CoulombFriction::CoulombFriction(char *matName,int matID) : ContactLaw(matName,matID)
{
	frictionCoeff = 0.0;			// <0 is stick
	frictionCoeffStatic = -1.;		// ignored if negative or if frictionCoeff < 0
	displacementOnly = 0.;			// >0 means displacementOnly, <=0 means tensile stress < abs(displacementOnly)
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
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&displacementOnly,gScaling,1.e6);
	}
	
	else if(strcmp(xName,"Dc")==0)
	{	input=DOUBLE_NUM;
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
	
	if(Dc<0.)
		cout << "   Stress found by perfect interface methods" << endl;
	else
	{	cout << "   Stress found by linear imperfect interface with Dc = ";
		cout << Dc*UnitsController::Scaling(1.e-6) << " " << UnitsController::Label(INTERFACEPARAM_UNITS) << endl;
	}
	if(displacementOnly>0.1)
		cout << "   Detection by only negative separation" << endl;
	else if(displacementOnly<0.)
	{	const char *label = UnitsController::Label(PRESSURE_UNITS);
		cout << "   Detection by negative separation and stress < " <<
			-displacementOnly*UnitsController::Scaling(1.e-6) << " " << label << endl;
	}
	else
		cout << "   Detection by negative separation and stress < 0" << endl;
}

#pragma mark CoulombFriction:Step Methods

// Change input momentum for stick (in delPi) to reflect frictional sliding contact law
//		(this called both my material contact and crack surface contact
// Input parameters (besides delPi)
//		norm = normal vector from material i to j
//		dotn = delPi.norm (precalculated)
//		deltaDotn = initial normal opening cod precalculated (=delta.norm)
//		mred = reduced mass
//		getHeating = true to calculated frictional heating term
//		contactArea = contact area, which is only needed by some laws
//		deltime = time step
//		delFi = force changed needed in post update calculations (only non-NULL in UPDATE_MOMENTUM_CALL)
//				and needed when frictional heating is activated
// Output
//		delPi change to reflect contact law
//		true is returned or false if decide now not in contact
//		*mredDelWf set to heat energy (actually mred*heat energy) (only if getHeating is true)
bool CoulombFriction::GetFrictionalDeltaMomentum(Vector *delPi,Vector *norm,double dotn,double deltaDotn,
							double *mredDelWf,double mred,bool getHeating,double contactArea,
							double deltime,Vector *delFi,NodalPoint *ndptr) const
{
	// indicate no frictional heat yet
	*mredDelWf=-1.;
	
	// stick and frictionless are easy and no heating
	if(frictionStyle==STICK)
	{	// stick conditions no change to delPi
		return true;
	}

	else if(frictionStyle==FRICTIONLESS)
	{	// Check contact at start or end of the interval
		
		// Scheme I - JANOSU-13-71 (scheme II OK, but scheme III terrible)
		double delEnd = deltaDotn;
		double fnaDt = dotn;
		if(delFi!=NULL)
			GetSeparationAndForce(delEnd,fnaDt,DotVectors(delFi,norm),deltime,mred,contactArea);
		
		// contact requires negative separation
		if(delEnd<0.)
		{	// Get stress cutoff
			double fnaDtMax = 0.;
			if(displacementOnly>0.1)
				fnaDtMax = fnaDt+1.;
			else if(displacementOnly<0.)
				fnaDtMax = -displacementOnly*contactArea*deltime;
			if(fnaDt<fnaDtMax)
			{	// frictionless contact - return normal component
				CopyScaleVector(delPi,norm,fnaDt);
				return true;
			}
		}

		// separated so no contact
		return false;
	}

	// Rest implements friction sliding
	// The initial delPi = (-N A dt) norm + (Sstick A dt) tang = dotn norm + dott tang
	// where N is normal traction (positive in compression), A is contact area, and dt is timestep
	
	// First verify contact during the time step and find the normal force
	double delEnd = deltaDotn;
	double fnaDt = dotn;
	if(delFi!=NULL)
		GetSeparationAndForce(delEnd,fnaDt,DotVectors(delFi,norm),deltime,mred,contactArea);
	
	// provisional setting for in contact
	bool inContact = false;
	if(delEnd<0.)
	{	// check stress or use just this displacement
		double fnaDtMax = 0.;
		if(displacementOnly>0.1)
			fnaDtMax = fnaDt+1.;
		else if(displacementOnly<0.)
			fnaDtMax = -displacementOnly*contactArea*deltime;
		if(fnaDt<fnaDtMax)
			inContact = true;
	}
	
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
	
	// Get frictional sliding force be Sslide Ac dt = f(N) Ac dt where NAcDt = -fnaDt
	double SslideAcDt = GetSslideAcDt(-fnaDt,dott,mred,contactArea,inContact,deltime);
	
	// if not in contact and did not find adhesive sticking, then done and no contact changes
	if(!inContact) return false;
	
	// if dott > Sslide Ac dt (which means Sstick>Sslide), then sliding, otherwise stick
	// but Sslide Ac dt <=0 means effectively frictionless
	if(SslideAcDt<=0.)
	{	// Normal stick condition
		CopyScaleVector(delPi,norm,fnaDt);
		//ZeroVector(delPi);		// fixes FricionL2 when displacement only because reverting to stress&displacement
	}
	else if(dott > SslideAcDt)
	{	// Normal stick condition
		CopyScaleVector(delPi,norm,fnaDt);

		// frictional terms only added in update momentum call
		if(delFi!=NULL)
			AddScaledVector(delPi,&tang,SslideAcDt);
		
		// get frictional heating term as friction work times reduced mass
		// As heat source need Del Wf/sec or divided by timestep*reduced mass
		// Note: only add frictional heating during momentum update (when friction
		//   force is appropriate) and only if frictional contact heat is enabled.
		if(getHeating)
		{	// delFi must be provided when getHeating is true (acceleration is delFi/mred, which is applied after return)
			// Vs alone is first order heating (i.e., *mredDelWf = Vs is first order method)
			double Vs = SslideAcDt*(dott-SslideAcDt);
			AddScaledVector(delFi, delPi, -1./deltime);
			double AsDt = SslideAcDt*DotVectors(delFi,&tang)*deltime;
			if(AsDt>Vs)
				*mredDelWf = 0.5*Vs*Vs/AsDt;
			else
				*mredDelWf = Vs - 0.5*AsDt;
			//*mredDelWf = Vs;							// revert to first order heating
		}
	}
	else
	{	// frictional stick - leave delPi as is
	}
	
	// still in contact
	return true;
}

// Get delta(dt) and interfacial force, which is only changed in post-momentum update contact (when delFi!=NULL)
// On input delEnd = delta(0) and fnaDt = Delta P_a.n = dn'
void CoulombFriction::GetSeparationAndForce(double &delEnd,double &fnaDt,double fn,double deltime,double mred,double Ac) const
{
	// perfect interface conditions (adjust separation, leave fnaDt alone)
	if(Dc<0.)
	{	delEnd += deltime*(fnaDt - 0.5*fn*deltime)/mred;
		return;
	}
	
	double m = mred/deltime;
	double d = Dc*Ac*deltime;
	
	// if frequency too high, use perfect interface limit
	if(d > 2.467401100272340*m)
	{	// treat as perfect  (adjust separation, leave fnaDt alone)
		delEnd += deltime*(fnaDt - 0.5*fn*deltime)/mred;
		return;
	}
	
	// linear imperfect interface analysis
	double sineTerm,sincosTerm;
	double phi = GetTerms(d,m,sineTerm,sincosTerm);
	
	double delZero = delEnd;
	double dnp = fnaDt;
	double fnDt = fn*deltime;
	
	fnaDt = 2.*(m*delEnd*(1-cos(phi)) + dnp*sineTerm) + fnDt*(2.*sincosTerm - 1.);
	delEnd = delZero*cos(phi) + (dnp*(1.-sineTerm) - fnDt*sincosTerm)/m;
}

// Return Sslide Ac dt = f(N) Ac dt
// Input is N Ac dt (and is always positive when in contact)
// If needed in the friction law, Ac is the contact area (not yet provided)
// The relative sliding speed after correcting the momentum will be (SStickAcDt-SslideAcDt)/mred
double CoulombFriction::GetSslideAcDt(double NAcDt,double SStickAcDt,double mred,
									  double contactArea,bool &inContact,double deltime) const
{
	// if not in contact return or if tension return zero friction force
	if(!inContact || NAcDt<0.) return 0.;
	
	// check static coefficient
	if(frictionCoeffStatic>0.)
	{	// If force to stick is less that static coefficient, then it will stick
		double Sstatic = frictionCoeffStatic*NAcDt;
		
		// returns force higher than static so contact law will pick static force
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

// True to ignore contact and revert to single velocity field
bool CoulombFriction::IgnoreContact(void) const { return false; }

// All interfaces need the law, if friction law needs it, must override and return true
bool CoulombFriction::ContactLawNeedsContactArea(void) const
{ 	return Dc>=0. || displacementOnly<0.;
}

// Return true is frictionless contact and no adhesion
bool CoulombFriction::IsFrictionless(void) const { return frictionStyle==FRICTIONLESS; }

// Return true is stick contact
bool CoulombFriction::IsStick(void) const { return frictionStyle==STICK; }

