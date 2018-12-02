/********************************************************************************
	ContactLaw.cpp
	nairn-mpm-fea

	Base class for all contact laws and itselt if subcalss of material
	By itself handle ignore contact

	Created by John Nairn, Oct 24, 2015.
	Copyright (c) 2015 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Materials/ContactLaw.hpp"
#include "Read_XML/MaterialController.hpp"
#include "Materials/CoulombFriction.hpp"
#include "Materials/LinearInterface.hpp"
#include "Nodes/CrackVelocityField.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark ContactLaw::Constructors and Destructors

// Constructor
ContactLaw::ContactLaw(char *matName,int matID) : MaterialBase(matName,matID)
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

#pragma mark ContactLaw::Initialization

// Read material properties by name (in xName). Set input to variable type
// (DOUBLE_NUM or INT_NUM) and return pointer to the class variable
// (cast as a char *)
// It gets called as a material, but really setting contact properties
char *ContactLaw::InputMaterialProperty(char *xName,int &input,double &gScaling)
{	return InputContactProperty(xName,input,gScaling);
}
char *ContactLaw::InputContactProperty(char *xName,int &input,double &gScaling)
{
	// This has no material properties, but ignore rho if there
	// ignore rho if there
	if(strcmp(xName,"rho")==0)
	{	input = NO_INPUT;
		return (char *)&rho;
	}
	
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
bool ContactLaw::GetFrictionalDeltaMomentum(Vector *delPi,Vector *norm,double dotn,double deltaDotn,double *mredDelWf,double mred,
											bool getHeating,double contactArea,double deltime,Vector *delFi) const
{
	return false;
}

// Contact handled here only for perfect interface parts (Dt or Dn < 0) by changing delPi
//		if done return false (but no used now?)
// Imperfect interfaces are calculated and always check if force is too high as determined by whether or not
//      the material's position is forced to pass the center of mass position by the calculated force
// Outputs are fImp, rawEnergy, and possible changed depPi. Rest are inputs, tandDel, deln, delt refer to cod vector
void ContactLaw::GetInterfaceForces(Vector *norm,Vector *fImp,double *rawEnergy,double surfaceArea,Vector *delPi,
										   double dotn,double mred,Vector *tangDel,double deln,double delt,double hperp) const
{}

// Return Sslide Ac dt = f(N) Ac dt
// For use by frictional contact only
// contactArea only provided if law reports it needs it in ContactLawNeedsContactArea()
// inContact will be true unless the law request to here about non-contact conditions
double ContactLaw::GetSslideAcDt(double NAcDt,double SStickAcDt,double mred,
								 double contactArea,bool &inContact,double deltime) const
{
	// law can change whether or not in contact, but most laws exit here if not in contact
	if(!inContact) return 0;
	
	// in contact, but no frictional force in contact law
	return 0.;
}

// Used by imperfect interfaces and maybe contact laws
// Get sineTerm = 1-sin phi/phi and sincosTerm = sin phi/phi -  (1-cos phi)/phi^2 stable even for phi near zero
// return phi
double ContactLaw::GetTerms(double d,double m,double &sineTerm,double &sincosTerm) const
{
	double phi = sqrt(d/m);
	if(phi<0.02)
	{	double phi2 = phi*phi;
		sineTerm = phi2/6.;				// = 1-sin phi/phi within 1e-10
		sincosTerm = 0.5 - 0.125*phi2;	// = sin phi/phi - (1-cos phi)/phi^2) within 1e-10
	}
	else
	{	double sineRatio = sin(phi)/phi;
		sineTerm = 1. - sineRatio;
		sincosTerm = sineRatio - (1.-cos(phi))/(phi*phi);
	}
	return phi;
}

#pragma mark ContactLaw::Accessors

// return material style
int ContactLaw::MaterialStyle(void) const { return CONTACT_MAT; }

// return unique, short name for this material
const char *ContactLaw::MaterialType(void) const { return "Contact Law Material"; }

// Calculate maximum wave speed for material in mm/sec.
double ContactLaw::WaveSpeed(bool threeD,MPMBase *mptr) const { return 1.e-12; }

// Is frictional contact and not ignoring contact
bool ContactLaw::IsFrictionalContact(void) const
{	if(IsImperfectInterface()) return false;
	return !IgnoreContact();
}

// All interfaces need the law, if friction law needs it, must override and return true
bool ContactLaw::ContactLawNeedsContactArea(void) const { return IsImperfectInterface(); }

// only this base law, rest should override and return false
bool ContactLaw::IgnoreContact(void) const { return true; }

// Describe type of imperfect interface
bool ContactLaw::IsImperfectInterface(void) const { return false; }

// Return true is frictionless contact and no adhesion
bool ContactLaw::IsFrictionless(void) const { return false; }

// Return true is stick contact
bool ContactLaw::IsStick(void) const { return false; }

#pragma mark ContactLaw::Class Methods

// Given old style input for friction,Dn, Dnc, and Dt, create appropriate contact law,
// add to contact law cache on material controller, and return a temporary ID
// if contactProps==NULL, create friction law with useFriction (cannot be for interface)
// throws std::bad_alloc
int ContactLaw::ConvertOldStyleToContactLaw(MaterialController *matCtrl,ContactInput *contactProps,double useFriction,const char *why)
{	// reading old style input
	
	// create contact law to recreate old style parameters
	ContactLaw *newContactLaw = NULL;
	char tempName[100];
	int numLaw = matCtrl->NumAutoContactLaws()+1;
	double newFriction = contactProps!=NULL ? contactProps->friction : useFriction;
	if(newFriction<-10.)
	{	sprintf(tempName,"Ignore Contact (Auto %d for %s)",numLaw,why);
		newContactLaw = new ContactLaw(tempName,CONTACTLAW);
	}
	else if(newFriction<0.)
	{	sprintf(tempName,"Stick Contact (Auto %d for %s)",numLaw,why);
		newContactLaw = new CoulombFriction(tempName,COULOMBFRICTIONLAW);
		((CoulombFriction *)newContactLaw)->SetFrictionCoeff(-1.);
	}
	else if(newFriction<10.)
	{	sprintf(tempName,"Frictional Contact (Auto %d for %s)",numLaw,why);
		newContactLaw = new CoulombFriction(tempName,COULOMBFRICTIONLAW);
		((CoulombFriction *)newContactLaw)->SetFrictionCoeff(newFriction);
	}
	else
	{	// Better have contactProps here
		sprintf(tempName,"Imperfect Interface (Auto %d for %s)",numLaw,why);
		newContactLaw = new LinearInterface(tempName,LINEARINTERFACELAW);
		if(contactProps!=NULL)
		{	if(contactProps->Dnc<-100.e6)
				contactProps->Dnc=contactProps->Dn;
			((LinearInterface *)newContactLaw)->SetParameters(contactProps->Dn,
															  contactProps->Dnc,contactProps->Dt);
		}
	}
	
	// cache to be new material and set finalLawID
	int finalLawID = 1000 + matCtrl->AddAutoContactLaw(newContactLaw);
	newContactLaw->autoID = finalLawID;
	return finalLawID;
}




