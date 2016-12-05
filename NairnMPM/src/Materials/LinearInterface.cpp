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
#include "NairnMPM_Class/MeshInfo.hpp"

extern double timestep;

#pragma mark LinearInterface::Constructors and Destructors

// Constructors
LinearInterface::LinearInterface() {}

LinearInterface::LinearInterface(char *matName) : ContactLaw(matName)
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
char *LinearInterface::InputMaterialProperty(char *xName,int &input,double &gScaling)
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
    return ContactLaw::InputMaterialProperty(xName,input,gScaling);
}

// Verify input properties do calculations; if problem return string with an error message
// Don't pass on to material base
const char *LinearInterface::VerifyAndLoadProperties(int np)
{
	// set Dnc if needed
	if(!hasSetDnc) Dnc = Dnt;
	
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
//		if done return false (but no used now?)
// Imperfect interfaces are calculated and always check if force is too high as determined by whether or not
//      the material's position is forced to pass the center of mass position by the calculated force
// Outputs are fImp, rawEnergy, and possible changed depPi. Rest are inputs, tandDel, deln, delt refer to cod vector
bool LinearInterface::GetInterfaceForcesForNode(Vector *norm,Vector *fImp,double *rawEnergy,
								double surfaceArea,Vector *delPi,double dotn,bool inContact,bool postUpdate,double mred,
								Vector *tangDel,double deln,double delt) const
{
	// initialize interfacial forces
    double trn=0.,trt=0.;
	
    // Convert delPi to shear momentum only (delPi - dotn (n) = dott (t)), then
	//   if shear force limited leave alone, otherwise set to zero
	//   if normal force limited, add normal back, otherwise leave alone
    AddScaledVector(delPi,norm,-dotn);
    
    // get shear (which is will always be positive because sign of tangent changes to accomodate it)
	// in g/(mm sec^2)
    if(Dt>=0.)
	{	// get force and compare to maximum force
		trt = Dt*delt*surfaceArea;
		double dott=sqrt(DotVectors(delPi,delPi));
		double maxFt = 2.*dott/timestep + 2.*mred*delt/(timestep*timestep);
        if(trt > maxFt)
		{	// limit force and retain delPi (which is for shear now)
            trt = 0.;
		}
        else
		{	// force OK, so remove shear momentum change now
            ZeroVector(delPi);
		}
    }
	// no 'else' needed because trt is zero, which means shear direction is flagged as perfect
	//     and retain delPi (which is already for shear stick)
    
	// get normal traction in g/(mm sec^2) - but different separated or in contact
	double maxFn;
	if(deln>0.)
	{	if(Dnt>=0.)
		{   // normal direction in tension is imperfect
			trn = Dnt*deln*surfaceArea;
			maxFn = 2.*dotn/timestep + 2.*mred*deln/(timestep*timestep);
			if(trn > maxFn)
			{	// limit force and add normal momentum change back into delPi
				trn = 0;
				AddScaledVector(delPi,norm,dotn);
			}
		}
		else
		{   // normal direction in tension flagged perfect
			// limit forces (trn already 0) and add normal momentum change back into delPi
			AddScaledVector(delPi,norm,dotn);
		}
	}
	else if(Dnc>=0.)
	{	// normal direction in compression is imperfect
		trn = Dnc*deln*surfaceArea;
		maxFn = 2.*dotn/timestep + 2.*mred*deln/(timestep*timestep);
		if(trn < maxFn)
		{	// limit negative force and add normal momentum change back into delPi
			trn = 0;
			AddScaledVector(delPi,norm,dotn);
		}	
	}
    else
    {   // normal direction is compression flagged as perfect
		// limit forces (trn already 0) and add normal momentum change back into delPi
        AddScaledVector(delPi,norm,dotn);
    }
    
    // find (trn n + trt t)*Ai for force in cartesian coordinates
    CopyScaleVector(fImp, norm, trn);
    AddScaledVector(fImp, tangDel, trt);
    
	// linear elastic energy
    *rawEnergy = 0.5*(trn*deln + trt*delt);
	
	return true;
}

// Get interface force for crack with internal interface
// da and db are displacements above and below the crack, and norm is normal vector (normalized)
// surfaceArea is contact surface area
// Output is force is fImp and renerge in rawEnergy
bool LinearInterface::GetCrackInterfaceForce(Vector *da,Vector *db,Vector *norm,double surfaceArea,double hperp,
												Vector *fImp,double *rawEnergy) const
{
	double dn,dt,trn = 0.,trt = 0.;
	
	if(Dnt>=0. || Dnc>=0.)
	{	// normal displacement
		dn = (db->x-da->x)*norm->x + (db->y-da->y)*norm->y;
		if(!mpmgrid.GetContactByDisplacements())
		{	// for efficiency used calculated dist
			dn -= mpmgrid.positionCutoff*hperp;
		}
		
		// Normal traction in g/(mm sec^2) - but different separated or in contact
		if(dn>0.)
		{	// normal direction in tension
			if(Dnt>=0.)
				trn = Dnt*dn*surfaceArea;
			else
			{	// interface perfect in tension, if also perfect in shear can exit
				if(Dt<0.) return false;
				dn = 0.;
			}
		}
		else
		{	// normal direction in compression
			if(Dnc>=0.)
				trn = Dnc*dn*surfaceArea;
			else
			{	// interface perfect in compression, if also perfect in shear can exit
				if(Dt<0.) return false;
				dn = 0.;
			}
		}
	}
	else
	{	// perfect in normal direction
		dn = 0.;
	}
	
	if(Dt>=0.)
	{	// transverse force
		dt = (db->x-da->x)*norm->y - (db->y-da->y)*norm->x;
		
		// transverse traction in g/(mm sec^2)
		trt = Dt*dt*surfaceArea;
	}
	else
	{	// perfect in normal direction
		dt = 0.;
	}
	
	// find trn n + trt t and finally normalize
	fImp->x = trn*norm->x + trt*norm->y;
	fImp->y = trn*norm->y - trt*norm->x;
	
	// total energy (not increment) is (1/2)(trn dnunnorm + trt dtunnorm)/(norm2*norm2) in g/sec^2
	// Use norm2 because |norm| for trn and trt and |norm| for dnunnorm and dtunnorm
	// units wiill be g/sec^2
	*rawEnergy = (trn*dn + trt*dt)/2.;
	
	return true;
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
bool LinearInterface::IsPerfectTangentialInterface(void) const { return Dt<0.; }
bool LinearInterface::IsPerfectNormalInterface(bool inContact) const
{	if( (!inContact && Dnt<0.) || (inContact && Dnc<0.) )
		return true;
	return false;
}
bool LinearInterface::IsPerfectNormalInterface(void) const { return Dnt<0. && Dnc>0.; }



