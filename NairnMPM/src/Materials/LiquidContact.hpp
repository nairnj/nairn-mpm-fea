/********************************************************************************
	LiquidContact.hpp
	nairn-mpm-fea
 
	Created by John Nairn, Feb 9, 2017.
	Copyright (c) 2017 John A. Nairn, All rights reserved.
 
	Dependencies
		CoulombFriction,hpp, ContactLaw.hpp, MaterialBase.hpp
 ********************************************************************************/

#ifndef LIQUIDCONTACT

#define LIQUIDCONTACT 64

#include "Materials/CoulombFriction.hpp"

class MaterialController;

class LiquidContact : public CoulombFriction
{
	public:
	
		// constructors and destructors
		LiquidContact(char *,int);
	
		// initialize
		virtual char *InputContactProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintContactLaw(void) const;
	
		// methods
		virtual double GetSslideAcDt(double,double,double,double,bool &,double) const;
	
		// accessors
		virtual const char *MaterialType(void) const;
		virtual bool ContactLawNeedsContactArea(void) const;
	
	protected:
		int liquidPhaseID;
		MaterialBase *liquidPhase;
};

#endif

