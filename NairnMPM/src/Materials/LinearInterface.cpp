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
		Dnc = Dnt;
	else if(DbleEqual(Dnc,Dnt))
		hasSetDnc = false;
	
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
												double dPDotn,double mred,Vector *tangDel,double deltaDotn,double deltaDott,double hperp) const
{
	// initialize interfacial forces and energy
    double trn=0.,trt=0.;
	*rawEnergy = 0.;
	double m = mred/timestep;
	double sineTerm,sincosTerm;
	
    // Convert delPi to shear momentum only (delPi - dPDotn (n) = dPDott (t)), then
	//   if shear force limited leave alone, otherwise set to zero
	//   if normal force limited, add normal back, otherwise leave alone
    AddScaledVector(delPi,norm,-dPDotn);

    // Get shear (which is will always be positive because sign of tangent changes to accomodate it)
	// Legacy units for force microN and displacement mm
    if(Dt>=0.)
	{	// get force and compare to maximum force

		// The following code is documented in contactetc.tex
		double dPDott=sqrt(DotVectors(delPi,delPi));
		double d = Dt*surfaceArea*timestep;
		
		// we limit to phi = sqrt(d/m) to less than pi/2 (pi^2/4 = 2.467401100272340)
		// May prefer to limit even futher, but never higher
		if(d > 2.467401100272340*m)
		{	// acceleration looks too high for stability, so revert to stick
			// set trt=0. (which initialized above) and keep delPi at stick conditions
		}
		else
		{	// exact solution method
			double phi = GetTerms(d,m,sineTerm,sincosTerm);
			double delFt = DotVectors(fImp,tangDel);

			// get force and final discontinuity
			trt = (2./timestep)*(m*deltaDott*(1-cos(phi)) + dPDott*sineTerm) + delFt*(2.*sincosTerm - 1.);
			double dut = deltaDott*cos(phi) + dPDott*(1.-sineTerm)/m - delFt*timestep*sincosTerm/m;

			// get total energy, because incremental approach does not seem to work
			*rawEnergy = 0.5*Dt*surfaceArea*dut*dut;
			
			// interface handled, so remove transverse stick condition
			ZeroVector(delPi);
		}
    }
	
	// Get normal traction - but different separated or in contact
	// Legacy units for force microN and displacement mm
	
	// if current perfect, stays perfect
	if((Dnt<0. && deltaDotn>=0.) || (Dnc<0. && deltaDotn<0.))
	{	// perfect interface
		// Set trn=0 (which is initialized above), and add normal momentum change back into delPi to make it perfect
		AddScaledVector(delPi,norm,dPDotn);
	}
	
	else
	{	// extract initial properties
		double Di, Df;
		if(deltaDotn>=0.)
		{	// initially opened
			Di = Dnt;
			Df = Dnc;
		}
		else
		{	Di = Dnc;
			Df = Dnt;
		}
		double d = Di*surfaceArea*timestep;
		
		// we limit to phi = sqrt(d/m) to less than pi/2 (pi^2/4 = 2.467401100272340)
		if(d > 2.467401100272340*m || Di<0.)
		{	// acceleration looks too high for stability, so revert to stick
			// Set trn=0 (which is initilized above),dd normal momentum change back into delPi to make it perfect
			AddScaledVector(delPi,norm,dPDotn);
		}
		else
		{	// exact solution method
			double phi = GetTerms(d,m,sineTerm,sincosTerm);
			double delFn = DotVectors(fImp,norm);
			
			// get force and final displacememt
			trn = (2./timestep)*(m*deltaDotn*(1-cos(phi)) + dPDotn*sineTerm) + delFn*(2.*sincosTerm - 1.);
			double dun = deltaDotn*cos(phi) + dPDotn*(1.-sineTerm)/m - delFn*timestep*sincosTerm/m;
			
			// done if no sign change
			if(dun*deltaDotn>=0. || !hasSetDnc)
			{	// get energy
				*rawEnergy += 0.5*Di*surfaceArea*dun*dun;
			}
			else
			{	// contact happens in the time step, we need to find the contact time
				double dn = dPDotn - delFn*timestep;
				double kap2 = Di*surfaceArea/mred;
				double kap = sqrt(kap2);
				
				// Find tc
				double A,B,C,tc,q,rootTerm;
				if(phi<1.e-5)
				{	// need special case for small phi
					A = delFn;
					B = 2.*dn;
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
				{	// These terms are in contactetc nores, but divided through by mred^2 kap^4
					double arg = deltaDotn - delFn/(mred*kap2);
					A = arg*arg + dn*dn/(mred*mred*kap2);
					B = 2.*delFn*dn/(mred*mred*kap*kap2);
					C = (delFn/(mred*kap2)+arg)*(delFn/(mred*kap2)-arg);
					
					// These comment forms are direct from notes
					//double arg = mred*kap2*deltaDotn - delFn;
					//A = arg*arg + kap2*dn*dn;
					//B = 2.*kap*delFn*dn;
					//C = (delFn+arg)*(delFn-arg);
					
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
				
				// find dn^(c) = dPDotn now to stick after contact
				// see contactetc results for first definition of dn^(c)
				phi = kap*tc;
				sineTerm = (phi<0.02) ? 1.-phi*phi/6 : sin(phi)/phi;
				dPDotn = dn*cos(phi) - mred*kap*deltaDotn*sin(phi) + delFn*tc*sineTerm;
				
				// post contact terms using Df
				kap2 = Df*surfaceArea/mred;
				kap = sqrt(kap2);
				
				// check if opposite direction is perfect
				if(Df<0. || kap2>2.467401100272340)
				{	// effective force to reach contact, but no energy because final discontinity is zero
					// found from contactetc notes with delta(delta t)=0
					trn = (2./timestep)*(m*deltaDotn + dPDotn) - delFn;
					
					// this version gets sum delPi + trn*timestep to equal final force, but poor results
					//trn = (2./timestep)*(m*deltaDotn - dPDotn) - delFn;
					
					// add momentum for post-contact stick
					//AddScaledVector(delPi,norm,dPDotn);
				}
				else
				{	// effective force for bilinear case
					// See contactetc notes - form using dn^(c) = dPDotn and not the expanded from
					phi = kap*timestep;
					if(phi<1.e-5)
					{	// for nearly debonded
						double dtTerm = 1. - tc/timestep;
						trn = (2./timestep)*(m*deltaDotn + dn - dPDotn*dtTerm) + delFn*(1 - dtTerm*dtTerm);
						
						// final discontinuity
						dtTerm = timestep - tc;
						dun = (dtTerm/mred)*(dPDotn + 0.5*delFn*dtTerm);
					}
					else
					{	// otherwise
						trn = (2./timestep)*(m*deltaDotn + dn - dPDotn*sin(kap*(timestep-tc))/phi)
									+ delFn*(1 - 2.*(1.-cos(kap*(timestep-tc)))/(phi*phi));
						
						// final discontinuity
						dun = (1./(mred*kap))*(dPDotn*sin(kap*(timestep-tc)) + delFn*(1-cos(kap*(timestep-tc)))/kap);
					}
					
					// energy from final dun
					*rawEnergy += 0.5*Df*surfaceArea*dun*dun;
				}
			}
		}
	}
	
    // find (trn n + trt t) for force in cartesian coordinates
    CopyScaleVector(fImp, norm, trn);
    AddScaledVector(fImp, tangDel, trt);
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

// True to ignore contact and revert to single velocity field
bool LinearInterface::IgnoreContact(void) const { return false; }

// True if model interface with tractions or false if handling contact
bool LinearInterface::IsImperfectInterface(void) const { return true; }



