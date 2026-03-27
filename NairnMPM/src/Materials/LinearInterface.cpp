/********************************************************************************
	LinearInterface.cpp
	nairn-mpm-fea

	Imperfect interface with linear traction law and not failure

	Created by John Nairn, Oct 24, 3015.
	Copyright (c) 2015 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/LinearInterface.hpp"
#include "System/UnitsController.hpp"

static int madeContact = 0;

extern double timestep;

#pragma mark LinearInterface::Constructors and Destructors

// Constructor
LinearInterface::LinearInterface(char *matName,int matID) : ContactLaw(matName,matID)
{
	Dnt = -UnitsController::Scaling(1.e6);
	Dnc = -UnitsController::Scaling(1.e6);
	Dt = -UnitsController::Scaling(1.e6);
	hasSetDnc = false;
}

#pragma mark LinearInterface::Initialization

// Read material properties by name (in xName). Set input to variable type
// (DOUBLE_NUM or INT_NUM) and return pointer to the class variable
// (cast as a char *)
char *LinearInterface::InputContactProperty(char *xName,int &input,double &gScaling)
{
    if(strcmp(xName,"Dn")==0 || strcmp(xName,"Dnt")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&Dnt,gScaling,1.e6);
	}
	
    else if(strcmp(xName,"Dnc")==0)
	{	input=DOUBLE_NUM;
		hasSetDnc = true;
		return UnitsController::ScaledPtr((char *)&Dnc,gScaling,1.e6);
	}
	
    else if(strcmp(xName,"Dt")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&Dt,gScaling,1.e6);
	}
	
	// does not all any from MaterialBase
    return ContactLaw::InputContactProperty(xName,input,gScaling);
}

// Verify input properties do calculations; if problem return string with an error message
// Don't pass on to material base
const char *LinearInterface::VerifyAndLoadProperties(int np)
{
	// set Dnc if needed or check if same
	if(!hasSetDnc)
	{	// set the same when linear in case used in code
		Dnc = Dnt;
	}
	else if(DbleEqual(Dnc,Dnt))
	{	// no lead to treat as bilinear
		hasSetDnc = false;
	}
	
	// must call super class
	return ContactLaw::VerifyAndLoadProperties(np);
}

// print contact law details to output window
void LinearInterface::PrintContactLaw(void) const
{
	const char *label = UnitsController::Label(INTERFACEPARAM_UNITS);
	PrintProperty("Dnt",Dnt*UnitsController::Scaling(1.e-6),label);
	PrintProperty("Dnc",Dnc*UnitsController::Scaling(1.e-6),label);
	PrintProperty("Dt",Dt*UnitsController::Scaling(1.e-6),label);
	cout << endl;
}

#pragma mark LinearInterface:Step Methods

// Contact handled here only for perfect interface parts (Dt or Dn < 0) by changing delPi
//		if done return false (i.e., if perfect interface with trn=trt=0
// Imperfect interfaces are calculated and always check if force is too high as determined by whether or not
//      the material's position is forced to pass the center of mass position by the calculated force
// On input, fImp = Delta f_a = ma Fc/Mc - Fa
// Inputs are: tangDel, deltaDotn, deltaDott and they refer to cod vector
// Outputs are fImp, rawEnergy, and possibly changed delPi.
void LinearInterface::GetInterfaceForces(Vector *norm,Vector *fImp,double *rawEnergy,double surfaceArea,Vector *delPi,
												double dPDotn,double mred,Vector *tangDel,double deltaDotn,double deltaDott,bool postUpdate) const
{
	// initialize interfacial forces and energy
    double trn=0.,trt=0.;
	*rawEnergy = 0.;
	double m = mred/timestep;
	double sineTerm,cosineTerm;
	double Cint2 = 2.25;
	
    // Convert delPi to shear momentum only (delPi - dPDotn (n) = dPDott (t)), then
	//   if shear force limited leave alone, otherwise set to zero
	//   if normal force limited, add normal back, otherwise leave alone
    AddScaledVector(delPi,norm,-dPDotn);

#pragma mark Shear Direction
	
    // Get shear (which is will always be positive because sign of tangent changes to accomodate it)
	// Legacy units for force microN and displacement mm
    if(Dt>=0.)
	{	// get force and compare to maximum force

		// The following code is documented in contactetc.tex
		double dPDott=sqrt(DotVectors(delPi,delPi));
		double d = Dt*surfaceArea*timestep;
		
		// Note that phi = sqrt(d*timestep/mred) = sqrt(d/m)
		// We limit phi to less than Cint or to d < Cint2*m
		// 2 is based on maximum time step for harmonic oscillator (lower might be better
		
		if(d > Cint2*m)
		{	// acceleration looks too high for stability, so revert to stick
			// leave trt=0. (initialized above) and keep delPi at stick conditions
		}
		else if(postUpdate)
		{	// exact solution method
			double phi = GetTerms(d,m,sineTerm,cosineTerm);
			double delFt = DotVectors(fImp,tangDel);

			// get force and final discontinuity
			double pt = dPDott - delFt*timestep;
			trt = Dt*surfaceArea*deltaDott*cosineTerm + 2.*pt*sineTerm + delFt*(1.-cosineTerm);
			double dut = deltaDott*cos(phi) + pt*(1.-sineTerm)/m + delFt*timestep*cosineTerm/(2.*m);

			// get total energy, because incremental approach does not seem to work
			*rawEnergy += 0.5*Dt*surfaceArea*dut*dut;
			
			// interface handled, so remove transverse stick condition
			ZeroVector(delPi);
		}
		else
		{	// interface will be handled in postUpdate
			// remove stick conditions now
			ZeroVector(delPi);
		}
    }

	// If handled shear than delPi = 0
	// If too stiff then delPi = dPDott (t)

#pragma mark Normal Direction
	// Get normal traction - but different separated or in contact
	// Legacy units for force microN and displacement mm
	
	// check for perfect directions
	double dnA = Dnt*surfaceArea*timestep;
	bool stiffN = Dnt<0. || dnA > Cint2*m ;
	double dnC = Dnc*surfaceArea*timestep;
	bool stiffC = Dnc<0. || dnC > Cint2*m ;

	// if both perfect then done, otherwise look for contact
	if(stiffN && stiffC)
	{	// perfect interface
		// Set trn=0 (which is initialized above), and add normal momentum change back into delPi to make it perfect
		AddScaledVector(delPi,norm,dPDotn);
	}
	
	else if(stiffC || stiffN)
	{	// one perfect direction, check the other (linear never here)
		double Di = stiffN ? Dnc : Dnt ;
		double delnSign = stiffN ? -1. : 1. ;
		double d = Di*surfaceArea*timestep;
			
		// exact solution method
		double phi = GetTerms(d,m,sineTerm,cosineTerm);
		double delFn = DotVectors(fImp,norm);

		// get force and final discontinuity
		double pn = dPDotn - delFn*timestep;
		trn = Di*surfaceArea*deltaDotn*cosineTerm + 2.*pn*sineTerm + delFn*(1.-cosineTerm);
		double dun = deltaDotn*cos(phi) + pn*(1.-sineTerm)/m + delFn*timestep*cosineTerm/(2.*m);
		
		// we are done if moves in nonstiff direction
		if(dun*delnSign>0. || dun==0.)
		{	// get total energy
			*rawEnergy += 0.5*Di*surfaceArea*dun*dun;
			
			// interface handled, leave normal removed
		}
		else
		{	// if moved wrong direction, use stiff direction
			// Set trn=0  and add normal momentum change back into delPi to make it perfect
			trn = 0.;
			AddScaledVector(delPi,norm,dPDotn);
		}
	}
	
	else
	{	// two non-stiff directions
		
		// extract initial properties
		double Di, Df, delnSign = 1.;
		if(deltaDotn>=0.)
		{	// initially opened
			Di = Dnt;
			Df = Dnc;
		}
		else
		{	// initially closed
			Di = Dnc;
			Df = Dnt;
		}
		
		// might need two tries if and only if deltaDotn==0.
		bool makesContact = false;
		double delFn,phi;
		while(true)
		{	double d = Di*surfaceArea*timestep;
			
			// exact solution method
			phi = GetTerms(d,m,sineTerm,cosineTerm);
			delFn = DotVectors(fImp,norm);

			// get force and final discontinuity
			double pn = dPDotn - delFn*timestep;
			trn = Di*surfaceArea*deltaDotn*cosineTerm + 2.*pn*sineTerm + delFn*(1.-cosineTerm);
			double dun = deltaDotn*cos(phi) + pn*(1.-sineTerm)/m + delFn*timestep*cosineTerm/(2.*m);
				
			// we are done if linear interface, start and stop positive or negative, if end at zero,
			// ... or if start at zero and ends with expected sign
			if(!hasSetDnc || dun*deltaDotn>0. || dun==0. || (deltaDotn==0. && dun*delnSign>=0.))
			{	// get total energy
				*rawEnergy += 0.5*Di*surfaceArea*dun*dun;
				
				// interface handled, leave normal removed
				break;
			}
			else if(deltaDotn==0. && delnSign>0.)
			{	// tension failed, try compression
				Di = Dnc;
				Df = Dnt;
				delnSign = -1.;
			}
			else
			{	// there is intersection
				makesContact = true;
				break;
			}
		}	// end two pass on initial loop

#pragma mark ...Handle Bilinear Step With Contact
		// If needed, handle time step with contact
		if(makesContact)
		{	// contact happens in this time step, we first need to find the contact time
			double pn = dPDotn - delFn*timestep;
			double kap2 = Di*surfaceArea/mred;
			double kap = sqrt(kap2);

#pragma mark ... Find Contact Time
			// Find tc
			double A,B,C,tc,q,rootTerm;
			if(phi<1.e-5)
			{	// need special case for small phi
				A = delFn;
				B = 2.*pn;
				C = 2.*mred*deltaDotn;
				rootTerm = B*B-4.*A*C;
				if(rootTerm<0.)
				{	// assume must be zero, but round off error
					q = -0.5*B;
				}
				else
					q= B>0 ? -0.5*(B+sqrt(rootTerm)) : -0.5*(B-sqrt(rootTerm));
				
				// we need root between 0 and timestep
				tc = q/A;
				if(tc<0. || tc>timestep) tc = C/q;
				if(tc<0. || tc>timestep) tc = 0.;
			}
			
			else
			{	// These terms are in contactetc notes, but divided through by mred^2 kap^4
				double arg = deltaDotn - delFn/(mred*kap2);
				A = arg*arg + pn*pn/(mred*mred*kap2);
				B = 2.*delFn*pn/(mred*mred*kap*kap2);
				double argF = delFn/(mred*kap2);
				C = argF*argF - arg*arg;
				
				// stable quadratic solution method
				rootTerm = B*B-4.*A*C;
				if(rootTerm<0.)
				{	// assume must be zero, but round off error (also never observed in calculations)
					q = -0.5*B;
				}
				else
					q= B>0 ? -0.5*(B+sqrt(rootTerm)) : -0.5*(B-sqrt(rootTerm));
				
				// we need root between 0 and 1
				double sinkt = q/A;
				if(sinkt<0. || sinkt>1.)
				{	// first invalid, try the second
					sinkt = C/q;
				}
				else
				{	double optSinkt = C/q;
					if(optSinkt>=0. && optSinkt<=1.)
					{	// two possible roots, used second if first too large
						if(sinkt > sin(phi))
							sinkt = optSinkt;
					}
				}
				
				if(sinkt<0. || sinkt>1.)
				{	// Never observed in calculations, but exit with one choice
					tc = 0.;
				}
				else
					tc = asin(sinkt)/kap;
			}
			if(madeContact==0)
			{	doutCritical("%FirstContact tc/dt Dinit",tc/timestep,Di);
#pragma omp atomic
				madeContact++;
			}

#pragma mark ... Post Contact Forces
			// find pn^(c) = pn shifted to contact time
			// see contactetc results for first definition of pn^(c)
			phi = kap*tc;
			sineTerm = (phi<0.02) ? 1.-phi*phi/6 : sin(phi)/phi;
			double pnc = pn*cos(phi) + (delFn - Di*surfaceArea*deltaDotn)*tc*sineTerm;
			
			// post contact terms using Df
			kap2 = Df*surfaceArea/mred;
			
			// check if opposite direction is perfect
			if(Df<0. || kap2>Cint2)
			{	// effective force to reach contact, but no energy because final discontinity is zero
				// found from contactetc notes with delta(delta t)=0
				//trn = (2./timestep)*(m*deltaDotn + dPDotn) - delFn;
				
				// post contact is stick, just impose on full time step
				trn = 0.;
				AddScaledVector(delPi,norm,dPDotn);
			}
			else if(postUpdate)
			{	// only need to finish up in post update step
				// effective force for bilinear case
				// See contactetc notes - form using dn^(c) = dPDotn and not the expanded from
				kap = sqrt(kap2);
				phi = kap*timestep;
				double dtTerm = timestep - tc;
				double dun;
				if(phi<1.e-5)
				{	// final discontinuity when kap zero (or close to it)
					// First phi terms are quadratic and ignored here
					// interface handled, so leave normal removed
					dun = (dtTerm/mred)*(dPDotn + 0.5*delFn*dtTerm);
				}
				else
				{	// final discontinuity
					dun = (1./(mred*kap))*(pnc*sin(kap*dtTerm) + delFn*(1-cos(kap*dtTerm))/kap);
				}
				
				// effective force
				trn = delFn + (2./timestep)*(m*(deltaDotn - dun) + pn);
				
				// energy from final dun
				*rawEnergy += 0.5*Df*surfaceArea*dun*dun;
				
				// interface handled, so leave normal removed
			}
			else
			{	// interface will be handled in postUpdate
				// set trn=0 (in case was tried) and leave normal removed
				trn = 0.;
			}
		}
	}
	
#pragma mark Final output
	
	// find (trn n + trt t) for force in cartesian coordinates
	CopyScaleVector(fImp, norm, trn);
	AddScaledVector(fImp, tangDel, trt);
	
	// delPi is zero or stick in either direction if needed
	
	// rawEnergy is summed if handled (or zero whan stick)
}

// Get sineTerm = 1-sin phi/phi and cosineTerm = 2*(1-cos phi)/phi^2 stable even for phi near zero
// return phi
double LinearInterface::GetTerms(double d,double m,double &sineTerm,double &cosineTerm) const
{
	double phi = sqrt(d/m);
	if(phi<0.02)
	{	double phi2 = phi*phi;
		sineTerm = phi2/6.;				// = 1-sin phi/phi within 1e-10
		cosineTerm = 1. - 0.5*sineTerm;	// = sin phi/phi - (1-cos phi)/phi^2) within 1e-10
	}
	else
	{	sineTerm = 1. - sin(phi)/phi;
		cosineTerm = 2.*(1.-cos(phi))/(phi*phi);
	}
	return phi;
}

#pragma mark LinearInterface::Accessors

// return unique, short name for this material
const char *LinearInterface::MaterialType(void) const { return "Linear Imperfect Interface"; }

// Set all three parameters
void LinearInterface::SetParameters(double newDn,double newDnc,double newDnt)
{	Dnt = newDn;
	Dnc = newDnc;
	Dt = newDnt;
	hasSetDnc = true;
}

// Inperfect interfaces handle contact
bool LinearInterface::IgnoreContact(void) const { return false; }

// True if model interface with tractions or false if handling contact
bool LinearInterface::IsImperfectInterface(void) const { return true; }



