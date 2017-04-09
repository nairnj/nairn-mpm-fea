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

// To get method from my 2013 paper, define PAPERIMPINTMETHOD below and also
// must define MANDMIMPINT in NodalPoint.hpp
//#define PAPERIMPINTMETHOD

extern double timestep;

#pragma mark NonlinearInterface::Constructors and Destructors

// Constructors
NonlinearInterface::NonlinearInterface() {}

NonlinearInterface::NonlinearInterface(char *matName) : LinearInterface(matName)
{
	order = 1;				// default to use higher-order method (0 is zeroth order, reset is 1)
}

#pragma mark LinearInterface::Initialization

// Read material properties by name (in xName). Set input to variable type
// (DOUBLE_NUM or INT_NUM) and return pointer to the class variable
// (cast as a char *)
char *NonlinearInterface::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
	if(strcmp(xName,"order")==0)
	{	input=INT_NUM;
		return (char *)(&order);
	}
	
	// does not all any from MaterialBase
	return LinearInterface::InputMaterialProperty(xName,input,gScaling);
}

// print contact law details to output window
void NonlinearInterface::PrintContactLaw(void) const
{
	LinearInterface::PrintContactLaw();
	
	// only options are zero and 0
	double dorder = order==0 ? 0. : 1. ;
	PrintProperty("Order",dorder,"");
	cout << endl;
}

#pragma mark LinearInterface:Step Methods

// Contact handled here only for perfect interface parts (Dt or Dn < 0) by changing delPi
//		if done return false (i.e., if perfect interface with trn=trt=0
// Imperfect interfaces are calculated and always check if force is too high as determined by whether or not
//      the material's position is forced to pass the center of mass position by the calculated force
// On input, fImp = Delta f_a = ma Fc/Mc - Fa
// Inputs are: tangDel, deltaDotn, deltaDott and they refer to cod vector
// Outputs are fImp, rawEnergy, and possibly changed depPi.
void NonlinearInterface::GetInterfaceForces(Vector *norm,Vector *fImp,double *rawEnergy,double surfaceArea,Vector *delPi,
										 double dPDotn,double mred,Vector *tangDel,double deltaDotn,double deltaDott,double hperp) const
{
	// initialize interfacial forces and energy
	double trn=0.,trt=0.;
	*rawEnergy = 0.;
	double m = mred/timestep;
	
	// Convert delPi to shear momentum only (delPi - dPDotn (n) = dPDott (t)), then
	//   if shear force limited leave alone, otherwise set to zero
	//   if normal force limited, add normal back, otherwise leave alone
	AddScaledVector(delPi,norm,-dPDotn);
	
	// Get shear (which is will always be positive because sign of tangent changes to accomodate it)
	// Legacy units for force microN and displacement mm
	
#ifdef PAPERIMPINTMETHOD
	
	if(Dt>0.)
	{	trt = Dt*deltaDott*surfaceArea;
		double dPDott=sqrt(DotVectors(delPi,delPi));
		double maxFt = 2.*(dPDott/timestep + mred*deltaDott/(timestep*timestep));
		if(trt > maxFt)
		{	// limit force and keep momentum change for tangential direction
			trt = 0;
			// keep delPi at stick conditions
		}
		else
		{	// all done, have force, and energy
			*rawEnergy += GetFtEnergy(deltaDott,trt);
			
			// remove shear stick
			ZeroVector(delPi);
		}
	}
	
#else
	
	// Check for numerical stability
	double d = GetFtPrime(deltaDott,surfaceArea)*timestep;
	int response = CheckDtStability(d,m);
	if(response!=STABLE)
	{	// acceleration looks too high for stability, so revert to stick
		trt = 0.;
		// keep delPi for stick conditions, discard for debond
		if(response==FORCE_DEBOND) ZeroVector(delPi);
	}
	else
	{	// zeroth order
		trt = GetFt(deltaDott,surfaceArea);

		if(order!=0)
		{	// new terms needed
			double dPDott=sqrt(DotVectors(delPi,delPi));
			double delFt = DotVectors(fImp,tangDel);
			
			// update discontinuity (to mean value over time step)
			double dut1 = deltaDott - trt*timestep/(6.*m) + 0.5*dPDott/m - delFt*timestep/(3.*m);
			
			// Recheck for numerical stability
			d = GetFtPrime(dut1,surfaceArea)*timestep;
			response = CheckDtStability(d,m);
			if(response!=STABLE)
			{	// acceleration looks too high for stability, so revert to stick
				trt = 0.;
				/// keep delPi for stick conditions, discard for debond
				if(response==FORCE_DEBOND) ZeroVector(delPi);
			}
			else
			{	// adjusted force
				trt = 0.5*(trt + GetFt(dut1,surfaceArea));
				
				// final discontinuity
				double dut = deltaDott - trt*timestep/(2.*m) + dPDott/m - delFt*timestep/(2.*m);
	
				// Recheck for numerical stability
				d = GetFtPrime(dut,surfaceArea)*timestep;
				response = CheckDtStability(d,m);
				if(response!=STABLE)
				{	// acceleration looks too high for stability, so revert to stick
					trt = 0.;
					// keep delPi for stick conditions, discard for debond
					if(response==FORCE_DEBOND) ZeroVector(delPi);
				}
				else
				{	// all done, have force, and energy
					*rawEnergy += GetFtEnergy(dut,trt);
					// remove shear stick
					ZeroVector(delPi);
				}
			}
		}
		else
		{	// all done, have force, and energy
			*rawEnergy += GetFtEnergy(deltaDott,trt);
			// remove shear stick
			ZeroVector(delPi);
		}
	}
#endif
	
	// Get normal traction

#ifdef PAPERIMPINTMETHOD
	
	// Get normal traction - but different separated or in contact
	// Legacy units for force microN and displacement mm
	if(deltaDotn>0.)
	{	if(Dnt>=0.)
		{   // normal direction in tension is imperfect
			trn = Dnt*deltaDotn*surfaceArea;
			double maxFn = 2.*(dPDotn/timestep + mred*deltaDotn/(timestep*timestep));
			if(trn > maxFn)
			{	// limit force and add normal momentum change back into delPi to make normal perfect
				trn = 0;
				AddScaledVector(delPi,norm,dPDotn);
			}
			else
				*rawEnergy += 0.5*trn*deltaDotn;
			
		}
		else
		{   // normal direction in tension flagged perfect
			// limit forces (trn already 0) and add normal momentum change back into delPi
			AddScaledVector(delPi,norm,dPDotn);
		}
	}
	else if(Dnc>=0.)
	{	// normal direction in compression is imperfect
		trn = Dnc*deltaDotn*surfaceArea;
		double maxFn = 2.*(dPDotn/timestep + mred*deltaDotn/(timestep*timestep));
		if(trn < maxFn)
		{	// limit negative force and add normal momentum change back into delPi to make normal perfect
			trn = 0;
			AddScaledVector(delPi,norm,dPDotn);
		}
		else
			*rawEnergy += 0.5*trn*deltaDotn;
	}
	else
	{   // normal direction is compression and flagged as perfect
		// limit forces (trn already 0) and add normal momentum change back into delPi
		AddScaledVector(delPi,norm,dPDotn);
	}
	
#else

	// Check for numerical stability
	d = GetFnPrime(deltaDotn,surfaceArea)*timestep;
	response = CheckDnStability(d,m,deltaDotn);
	if(response!=STABLE)
	{	// limit force and add normal momentum change back into delPi to make normal perfect
		trn = 0.;
		// add back for stick, keep zero for debond
		if(response==FORCE_STICK) AddScaledVector(delPi,norm,dPDotn);
	}
	else
	{	// zeroth order
		trn = GetFn(deltaDotn,surfaceArea);
		
		if(order!=0)
		{	// new terms needed
			double delFn = DotVectors(fImp,norm);
		
			// update discontinuity
			double dun1 = deltaDotn - trn*timestep/(6.*m) + 0.5*dPDotn/m - delFn*timestep/(3.*m);
			
			// Recheck for numerical stability
			d = GetFtPrime(dun1,surfaceArea)*timestep;
			response = CheckDnStability(d,m,deltaDotn);
			if(response!=STABLE)
			{	// limit force and add normal momentum change back into delPi to make normal perfect
				trn = 0.;
				// add back for stick, keep zero for debond
				if(response==FORCE_STICK) AddScaledVector(delPi,norm,dPDotn);
			}
			else
			{	// adjusted force
				trn = 0.5*(trn + GetFn(dun1,surfaceArea));
				
				// final discontinuity
				double dun = deltaDotn - trn*timestep/(2.*m) + dPDotn/m - delFn*timestep/(2.*m);
				
				// Recheck for numerical stability
				d = GetFtPrime(dun,surfaceArea)*timestep;
				response = CheckDnStability(d,m,deltaDotn);
				if(response!=STABLE)
				{	/// limit force and add normal momentum change back into delPi to make normal perfect
					trn = 0.;
					// add back for stick, keep zero for debond
					if(response==FORCE_STICK) AddScaledVector(delPi,norm,dPDotn);
				}
				else
				{	// all done, have force, and energy
					*rawEnergy += GetFnEnergy(dun,trn);
				}
			}
		}
		else
		{	// all done, have force, and energy
			*rawEnergy += GetFnEnergy(deltaDotn,trn);
		}
	}
#endif
	
	// find (trn n + trt t) for force in cartesian coordinates
	CopyScaleVector(fImp, norm, trn);
	AddScaledVector(fImp, tangDel, trt);
}

// Nonlinear subclasses only need to support follow methods for given interface law.
// One set for tangential and one set for normal

// check stability
int NonlinearInterface::CheckDtStability(double d,double m) const
{	// stick if set to perfect
	if(Dt<0.) return FORCE_STICK;
	
	// see if stable (less then pi^2/2)
	if(fabs(d) <= 2.467401100272340*m) return STABLE;
	
	// instable, but should it debond or stick
	return d>0. ? FORCE_STICK : FORCE_DEBOND ;
}

// tangential force
double NonlinearInterface::GetFt(double delta,double area) const
{	return Dt*area*delta;
}

// tangential stiffness
double NonlinearInterface::GetFtPrime(double delta,double area) const
{	return Dt*area;
}

// tangential energy
double NonlinearInterface::GetFtEnergy(double dut,double trt) const
{	return 0.5*trt*dut;
}

// check stability
int NonlinearInterface::CheckDnStability(double d,double m,double delta) const
{	// stick if set to perfect
	if(delta>=0.)
	{	if(Dnt<0.) return FORCE_STICK;
	}
	else if(Dnc<0.)
		return FORCE_STICK;
	
	// see if stable (less then pi^2/2)
	if(fabs(d) <= 2.467401100272340*m) return STABLE;
	
	// instable, but should it debond or stick
	return d>0. ? FORCE_STICK : FORCE_DEBOND ;
}

// normal force
double NonlinearInterface::GetFn(double delta,double area) const
{	return delta>0 ? Dnt*area*delta : Dnc*area*delta;
}

// tangential stiffness
double NonlinearInterface::GetFnPrime(double delta,double area) const
{	return delta>0 ? Dnt*area : Dnc*area;
}

// tangential energy
double NonlinearInterface::GetFnEnergy(double dun,double trn) const
{	return 0.5*trn*dun;
}

#pragma mark LinearInterface::Accessors

// return unique, short name for this material
const char *NonlinearInterface::MaterialType(void) const { return "Nonlinear Imperfect Interface"; }
