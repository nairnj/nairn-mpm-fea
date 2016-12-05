/********************************************************************************
	ContactLaw.cpp
	nairn-mpm-fea

	Base class for all contact laws and itselt if subcalss of material
	By itself handle ignore contact

	Created by John Nairn, Oct 24, 3015.
	Copyright (c) 2015 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/ContactLaw.hpp"
#include "Read_XML/MaterialController.hpp"
#include "Materials/CoulombFriction.hpp"
#include "Materials/LinearInterface.hpp"

#pragma mark NewMaterial::Constructors and Destructors

// Constructors
ContactLaw::ContactLaw() {}

ContactLaw::ContactLaw(char *matName) : MaterialBase(matName)
{
	// initialized and only used to support old style contact input
	contactProps.friction=0.;				// crack contact friction
	contactProps.Dn=-1.;					// prefect in tension by default
	contactProps.Dnc=-101.e6;				// <-100e6 means not set and should be set same as Dn
	contactProps.Dt=-1.;					// perfect in shear by default
	contactProps.otherMatID=-1;
	contactProps.contactLawID=-1;
	autoID = -1;
}

#pragma mark NewMaterial::Initialization

// Read material properties by name (in xName). Set input to variable type
// (DOUBLE_NUM or INT_NUM) and return pointer to the class variable
// (cast as a char *)
char *ContactLaw::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
	// This has no material properties
	
	// does not all any from MaterialBase
    return((char *)NULL);
}

// Verify input properties do calculations; if problem return string with an error message
// Don't pass on to material base
const char *ContactLaw::VerifyAndLoadProperties(int np)
{
	// check properties
	
	// must call super class
	return NULL;
}

// print to output window
void ContactLaw::PrintMaterial(int num) const
{
	// material name
    cout << "Material " << num << ": " << name << endl;
    cout << "     " << MaterialType() << " with:" << endl;
	PrintContactLaw();
}

// print contact law details to output window
void ContactLaw::PrintContactLaw(void) const
{
	cout << "Contact nodes revert to center of mass velocity field" << endl;
}

#pragma mark ContactLaw:Step Methods

// Adjust change in momentum for frictional contact in tangential direction
// If has component of tangential motion, calculate force depending on whether it is sticking or sliding
// When frictional sliding, find tangential force (times dt) and set flag, if not set flag false
// When friction heating is on, set Ftdt term and set hasFriction to true
//     (hasFriction (meaning has frictional heating value) must be initialized to false when called)
// contactArea only provided if frictional law needs it
// Normally in contact when called, but some laws might want call even when not in contact. If return
//		value is false, the contact should be treated as no contact
bool ContactLaw::GetFrictionalDeltaMomentum(Vector *delPi,Vector *norm,double dotn,double *mredDE,double mred,
											bool getHeating,double contactArea,bool inContact,double deltime,Vector *at) const
{
	return false;
}

// Contact handled here only for perfect interface parts (Dt or Dn < 0) by changing delPi
//		if done return false (but no used now?)
// Imperfect interfaces are calculated and always check if force is too high as determined by whether or not
//      the material's position is forced to pass the center of mass position by the calculated force
// Outputs are fImp, rawEnergy, and possible changed depPi. Rest are inputs, tandDel, deln, delt refer to cod vector
bool ContactLaw::GetInterfaceForcesForNode(Vector *norm,Vector *fImp,double *rawEnergy,
								double surfaceArea,Vector *delPi,double dotn,bool inContact,bool postUpdate,double mred,
								Vector *tangDel,double deln,double delt) const
{
	return false;
}

// Get interface force for crack with internal interface
// da and db are displacements above and below the crack, and norm is normal vector (normalized)
// surfaceArea is contact surface area
// Output is force is fImp and renerge in rawEnergy
bool ContactLaw::GetCrackInterfaceForce(Vector *da,Vector *db,Vector *norm,double surfaceArea,double dist,
											 Vector *fImp,double *rawEnergy) const
{
	return false;
}

// Return Sslide Ac dt = f(N) Ac dt
// For use by frictional contact only
// contactArea only provided is law reports it needs it in FrictionLawNeedsContactArea()
// inContact will be true unless the law request to here about non-contact conditions
double ContactLaw::GetSslideAcDt(double NAcDt,double SStickAcDt,double Ac,
								 double mred,double contactArea,bool &inContact,double deltime) const
{	return 0.;
}

#pragma mark ContactLaw::Accessors

// return material style
int ContactLaw::MaterialStyle(void) const { return CONTACT_MAT; }

// return unique, short name for this material
const char *ContactLaw::MaterialType(void) const { return "Contact Law Material"; }

// Calculate maximum wave speed for material in mm/sec.
double ContactLaw::WaveSpeed(bool threeD,MPMBase *mptr) const { return 1.e-12; }

// Are calculations done depending on if found to be inContact
bool ContactLaw::ContactIsDone(bool inContact) const
{	if(IsImperfectInterface()) return false;
	if(!inContact) return true;
	return false;
}

// Give details aobut frictional contact
bool ContactLaw::IsFrictionalContact(void) const { return !IsImperfectInterface(); }
bool ContactLaw::FrictionLawNeedsContactArea(void) const { return false; }
bool ContactLaw::IgnoreContact(void) const { return true; }

// Describee type of imperfect interface
bool ContactLaw::IsImperfectInterface(void) const { return false; }
bool ContactLaw::IsPerfectTangentialInterface(void) const { return false; }
bool ContactLaw::IsPerfectNormalInterface(bool inContact) const { return false; }
bool ContactLaw::IsPerfectNormalInterface(void) const { return false; }

#pragma mark ContactLaw::Class Methods

// Given old style input for friction,Dn, Dnc, and Dt, create appropriate contact law,
// add to contact law cache on material controller, and return a temporary ID
// if contactProps==NULL, create friction law with useFriction (cannot be for interface)
// throws std::bad_alloc
int ContactLaw::ConvertOldStyleToContactLaw(MaterialController *matCtrl,ContactInput *contactProps,double useFriction)
{	// reading old style input
	
	// create contact law to recreate old style parameters
	ContactLaw *newContactLaw = NULL;
	char tempName[100];
	int numLaw = matCtrl->NumAutoContactLaws()+1;
	double newFriction = contactProps!=NULL ? contactProps->friction : useFriction;
	if(newFriction<-10.)
	{	sprintf(tempName,"Ignore Contact (Auto %d)",numLaw);
		newContactLaw = new ContactLaw(tempName);
	}
	else if(newFriction<0.)
	{	sprintf(tempName,"Stick Contact (Auto %d)",numLaw);
		newContactLaw = new CoulombFriction(tempName);
		((CoulombFriction *)newContactLaw)->SetFrictionCoeff(-1.);
	}
	else if(newFriction<10.)
	{	sprintf(tempName,"Frictional Contact (Auto %d)",numLaw);
		newContactLaw = new CoulombFriction(tempName);
		((CoulombFriction *)newContactLaw)->SetFrictionCoeff(newFriction);
	}
	else
	{	// Better have contactProps here
		sprintf(tempName,"Imperfect Interface (Auto %d)",numLaw);
		newContactLaw = new LinearInterface(tempName);
		if(contactProps->Dnc<-100.e6)
			contactProps->Dnc=contactProps->Dn;
		((LinearInterface *)newContactLaw)->SetParameters(contactProps->Dn,
														contactProps->Dnc,contactProps->Dt);
	}
	
	// cache to be new material and set finalLawID
	int finalLawID = 1000 + matCtrl->AddAutoContactLaw(newContactLaw);
	newContactLaw->autoID = finalLawID;
	return finalLawID;
}




