/********************************************************************************
	NoninearInterface.cpp
	nairn-mpm-fea
 
	Created by John Nairn, 3/30/2017.
	Copyright (c) 2017 John A. Nairn, All rights reserved.
 
	The class is really an abstract class to serve as parent to all
	non-linear interface laws. It can be used (by #59 only) it does
	linear interface. A directive can make linear mode identical to
	old interface from my paper while LinearInterface class now use
	harmonic oscillator analysis.
********************************************************************************/

#include "stdafx.h"
#include "Materials/NonlinearInterface.hpp"
#include "System/UnitsController.hpp"

// hack to pick 0th, 1st, or 2nd order
#define NL_ORDER 2

extern double timestep;

#pragma mark NonlinearInterface::Constructors and Destructors

// Constructor
NonlinearInterface::NonlinearInterface(char *matName,int matID) : LinearInterface(matName,matID)
{
	stylen = NL_LINEAR_INTERFACE;
	stylet = NL_LINEAR_INTERFACE;
	peakn = 0.;
	peakt = 0.;
}

#pragma mark LinearInterface::Initialization

// Read material properties by name (in xName). Set input to variable type
// (DOUBLE_NUM or INT_NUM) and return pointer to the class variable
// (cast as a char *)
char *NonlinearInterface::InputContactProperty(char *xName,int &input,double &gScaling)
{
	
	if(strcmp(xName,"normal_shape")==0)
	{	input=INT_NUM;
		return (char *)(&stylen);
	}
	
	else if(strcmp(xName,"tangential_shape")==0)
	{	input=INT_NUM;
		return (char *)(&stylet);
	}

	else if(strcmp(xName,"Npeak")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&peakn,gScaling,1.e6);
	}
	
	else if(strcmp(xName,"Tpeak")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&peakt,gScaling,1.e6);
	}
	
	// does not all any from MaterialBase
	return LinearInterface::InputContactProperty(xName,input,gScaling);
}

// Verify input properties do calculations; if problem return string with an error message
// Don't pass on to material base
const char *NonlinearInterface::VerifyAndLoadProperties(int np)
{
	if(stylet<NL_LINEAR_INTERFACE || stylet>=MAX_STYLES)
		return "Invalid nonlinear shape provided to tangential direction";
	
	if(stylen<NL_LINEAR_INTERFACE || stylen>=MAX_STYLES)
		return "Invalid nonlinear shape provided to tangential direction";
	
	const char *parent = LinearInterface::VerifyAndLoadProperties(np);
	if(parent!=NULL) return parent;
	
	// reject bilinear
	if(stylen==NL_LINEAR_INTERFACE && hasSetDnc)
		return "Must use LinearInterface to model bilinear interface with Dnc!=Dnt";
	
	// get Morse properties
	if(stylet==MORSE_POTENTIAL)
	{	// require Dt and peakt>0}
		if(Dt<=0. || peakt<=0.)
			return "Morse potential in tangential direction requires Dt>0 and Tpeak>0";
		Det = 8.*peakt*peakt/Dt;
		alphat = Dt/(4.*peakt);
	}
	if(stylen==MORSE_POTENTIAL)
	{	// require Dnt and peakn>0}
		if(Dnt<=0. || peakn<=0.)
			return "Morse potential in normal direction requires Dn>0 and Npeak>0";
		Den = 8.*peakn*peakn/Dnt;
		alphan = Dnt/(4.*peakn);
	}

	return NULL;
}

// print contact law details to output window
void NonlinearInterface::PrintContactLaw(void) const
{
	const char *label = UnitsController::Label(INTERFACEPARAM_UNITS);
	
	// normal
	switch(stylen)
	{	case NL_LINEAR_INTERFACE:
			cout << "Normal direction is linear:" << endl;
			PrintProperty("Dn",Dnt*UnitsController::Scaling(1.e-6),label);
			cout << endl;
			break;
		case MORSE_POTENTIAL:
			cout << "Normal direction is Morse potential" << endl;
			PrintProperty("Dn",Dnt*UnitsController::Scaling(1.e-6),label);
			PrintProperty("Npeak",peakn*UnitsController::Scaling(1.e-6),UnitsController::Label(PRESSURE_UNITS));
			PrintProperty("[umax]",(log(2.)/alphan),UnitsController::Label(CULENGTH_UNITS));
			cout << endl;
			break;
	}

	// tangential
	switch(stylet)
	{	case NL_LINEAR_INTERFACE:
			cout << "Tangential direction is linear:" << Dt*UnitsController::Scaling(1.e-6);
			PrintProperty("Dt",Dnt*UnitsController::Scaling(1.e-6),label);
			cout << endl;
			break;
		case MORSE_POTENTIAL:
			cout << "Tangential direction is Morse potential" << endl;
			PrintProperty("Dt",Dt*UnitsController::Scaling(1.e-6),label);
			PrintProperty("Tpeak",peakn*UnitsController::Scaling(1.e-6),UnitsController::Label(PRESSURE_UNITS));
			PrintProperty("[umax]",(log(2.)/alphat),UnitsController::Label(CULENGTH_UNITS));
			cout << endl;
			break;
	}
}

#pragma mark NonlinearInterface:Step Methods

// Contact handled here only for perfect interface parts (Dt or Dn < 0) by changing delPi
//		if done return false (i.e., if perfect interface with trn=trt=0
// Imperfect interfaces are calculated and always check if force is too high as determined by whether or not
//      the material's position is forced to pass the center of mass position by the calculated force
// On input, fImp = Delta f_a = ma Fc/Mc - Fa
// Inputs are: tangDel, deltaDotn, deltaDott and they refer to cod vector
// Outputs are fImp, rawEnergy, and possibly changed delPi.
void NonlinearInterface::GetInterfaceForces(Vector *norm,Vector *fImp,double *rawEnergy,double surfaceArea,Vector *delPi,
										 double dPDotn,double mred,Vector *tangDel,double deltaDotn,double deltaDott,bool postUpdate) const
{
	// initialize interfacial forces and energy
	double trn=0.,trt=0.;
	*rawEnergy = 0.;
	double m = mred/timestep;
	
	// Convert delPi to shear momentum only (delPi - dPDotn (n) = dPDott (t)), then
	//   if shear force limited leave alone, otherwise set to zero
	//   if normal force limited, add normal back, otherwise leave alone
	AddScaledVector(delPi,norm,-dPDotn);

#pragma mark Shear Direction

	// Get shear (which is will always be positive because sign of tangent changes to accomodate it)
	// Legacy units for force microN and displacement mm
	
	// Check for numerical stability: (d/m) = F'dt^2/mred = phi^2 (see contactetc notes)
	double d = GetFtPrime(deltaDott,surfaceArea)*timestep;
	if(CheckDtStability(d,m)!=STABLE)
	{	// acceleration looks too high for stability, so revert to stick
		// leave trt=0. (initialized above) and keep delPi at stick conditions
	}
	else if(postUpdate)
	{	// second order
		double dPDott=sqrt(DotVectors(delPi,delPi));
		double delFt = DotVectors(fImp,tangDel);
		double phi2 = d/m;
		double Ft0 = GetFt(deltaDott,surfaceArea);
		double fdiff = delFt - Ft0;
		double pt = dPDott - delFt*timestep;
		
		// force and discontinuity
#if NL_ORDER==2
		trt = Ft0 + pt*phi2/(3*timestep) + fdiff*phi2/12.;
		double dut = deltaDott + (pt/m)*(1.-phi2/6.) + (timestep*fdiff/(2.*m))*(1.-phi2/12.);
#elif NL_ORDER==1
		trt = Ft0 + pt*phi2/(3*timestep);
		double dut = deltaDott + (pt/m)*(1.-phi2/6.) + (timestep*fdiff/(2.*m));
#else
		trt = Ft0;
		double dut = deltaDott + (pt/m) + (timestep*fdiff/(2.*m));
#endif
		
		// get total energy
		*rawEnergy += GetFtEnergy(dut,surfaceArea);
		
		// interface handled, so remove transverse stick condition
		ZeroVector(delPi);
	}
	else
	{	// interface will be handled in postUpdate
		// remove stick conditions now
		ZeroVector(delPi);
	}

	// If handled shear than delPi = 0
	// If too stiff then delPi = dPDott (t)

#pragma mark Normal Direction
	
	// Check for numerical stability: (d/m) = F'dt^2/mred = phi^2 (see contactetc notes)
	d = GetFnPrime(deltaDotn,surfaceArea)*timestep;
	if(CheckDnStability(d,m,deltaDotn)!=STABLE)
	{	// acceleration looks too high for stability, so revert to stick
		// Set trt=0 (which is initialized above), and add normal momentum change back into delPi to make it perfect
		AddScaledVector(delPi,norm,dPDotn);
	}
	else if(postUpdate)
	{	// second order
		double delFn = DotVectors(fImp,norm);
		double phi2 = d/m;
		double Fn0 = GetFn(deltaDotn,surfaceArea);
		double fdiff = delFn - Fn0;
		double pn = dPDotn - delFn*timestep;
		
		// force and discontinuity
#if NL_ORDER==2
		trn = Fn0 + pn*phi2/(3*timestep) + fdiff*phi2/12.;
		double dun = deltaDotn + (pn/m)*(1.-phi2/6.) + (timestep*fdiff/(2.*m))*(1.-phi2/12.);
#elif NL_ORDER==1
		trn = Fn0 + pn*phi2/(3*timestep);
		double dun = deltaDotn + (pn/m)*(1.-phi2/6.) + (timestep*fdiff/(2.*m));
#else
		trn = Fn0;
		double dun = deltaDotn + (pn/m) + (timestep*fdiff/(2.*m));
#endif

		// get total energy
		*rawEnergy += GetFnEnergy(dun,surfaceArea);
		
		// interface handled, leave normal stick removed
		// both in postUpdate and if not postUpdate
	}

#pragma mark Final output
	
	// find (trn n + trt t) for force in cartesian coordinates
	CopyScaleVector(fImp, norm, trn);
	AddScaledVector(fImp, tangDel, trt);
	
	// delPi is zero or stick in either direction if needed
	
	// rawEnergy is summed if handled (or zero whan stick)
}

#pragma mark NonlinearInterface::Law Implementation

// Nonlinear subclasses only need to support follow methods for given interface law.
// One set for tangential and one set for normal.

// check stability
int NonlinearInterface::CheckDtStability(double d,double m) const
{
	if(stylet==NL_LINEAR_INTERFACE)
	{	// stick if set to perfect
		if(Dt<0.) return FORCE_STICK;
	}
	
	// see if stable (d/m <= Cint^2)
	return d <= 2.25*m ? STABLE : FORCE_STICK ;
}

// tangential force
double NonlinearInterface::GetFt(double delta,double area) const
{
	double Ft = 0.;
	switch(stylet)
	{	case NL_LINEAR_INTERFACE:
			Ft = Dt*area*delta;
			break;
		case MORSE_POTENTIAL:
			Ft = MorseForce(delta,alphat,Det)*area;
			break;
		default:
			break;
	}
	return Ft;
}

// tangential stiffness
double NonlinearInterface::GetFtPrime(double delta,double area) const
{
	double Ftp = 0.;
	switch(stylet)
	{	case NL_LINEAR_INTERFACE:
			Ftp = Dt*area;
			break;
		case MORSE_POTENTIAL:
			Ftp = MorseSlope(delta,alphat,Dt)*area;
			break;
		default:
			break;
	}
	return Ftp;
}

// tangential energy
double NonlinearInterface::GetFtEnergy(double dut,double area) const
{
	double energy = 0.;
	switch(stylet)
	{	case NL_LINEAR_INTERFACE:
			energy = 0.5*Dt*area*dut*dut;
			break;
		case MORSE_POTENTIAL:
			energy = MorseEnergy(dut,alphat,Det)*area;
			break;
		default:
			break;
	}
	return energy;
}

// check stability
int NonlinearInterface::CheckDnStability(double d,double m,double delta) const
{
	if(stylen==NL_LINEAR_INTERFACE)
	{	// stick if set to perfect
		if(Dnt<0.) return FORCE_STICK;
	}
	
	// see if stable (d/m <= Cint^2)
	return d <= 2.25*m ? STABLE : FORCE_STICK ;
}

// normal force
double NonlinearInterface::GetFn(double delta,double area) const
{
	double Fn = 0.;
	switch(stylen)
	{	case NL_LINEAR_INTERFACE:
			Fn = Dnt*area*delta;
			break;
		case MORSE_POTENTIAL:
			Fn = MorseForce(delta,alphan,Den)*area;
			break;
		default:
			break;
	}
	return Fn ;
}

// tangential stiffness
double NonlinearInterface::GetFnPrime(double delta,double area) const
{
	double Fnp = 0.;
	switch(stylen)
	{	case NL_LINEAR_INTERFACE:
			Fnp = Dnt*area;
			break;
		case MORSE_POTENTIAL:
			Fnp = MorseSlope(delta,alphan,Dnt)*area;
			break;
		default:
			break;
	}
	return Fnp ;
}

// tangential energy
double NonlinearInterface::GetFnEnergy(double dun,double area) const
{
	double energy = 0.;
	switch(stylen)
	{	case NL_LINEAR_INTERFACE:
			energy = 0.5*Dnt*area*dun*dun;
			break;
		case MORSE_POTENTIAL:
			energy = MorseEnergy(dun,alphan,Den)*area;
			break;
		default:
			break;
	}
	return energy;
}

#pragma mark NonlinearInterface::MorsePotential

// Get force from Morse potential
double NonlinearInterface::MorseForce(double x,double alpha,double De) const
{
	double arg = exp(-alpha*x);
	return 2.*De*alpha*arg*(1.-arg);
}

// Get slope of force from Morse potential
double NonlinearInterface::MorseSlope(double x,double alpha,double D) const
{
	double arg = exp(-alpha*x);
	return D*arg*(2.*arg-1.);
}

// Get energy from Morse potential
double NonlinearInterface::MorseEnergy(double x,double alpha,double De) const
{
	double arg = 1. - exp(-alpha*x);
	return De*arg*arg;
}

#pragma mark NonlinearInterface::Accessors

// return unique, short name for this material
const char *NonlinearInterface::MaterialType(void) const { return "Nonlinear Imperfect Interface"; }
