/********************************************************************************
	CoulombFriction.hpp
	nairn-mpm-fea

	Created by John Nairn, Oct 24, 2015.
	Copyright (c) 2015 John A. Nairn, All rights reserved.

	Dependencies
		ContactLaw.hpp,MaterialBase.hpp
********************************************************************************/

#ifndef COULOMBFRICTIONLAW

#define COULOMBFRICTIONLAW 61

#include "Materials/ContactLaw.hpp"

class CoulombFriction : public ContactLaw
{
	public:
	
		// constructors and destructors
		CoulombFriction(char *,int);
	
		// initialize
		virtual char *InputContactProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintContactLaw(void) const;
	
		// methods
		virtual bool GetFrictionalDeltaMomentum(Vector *,Vector *,double,double,double *,double,bool,double,double,Vector *) const;
		virtual void GetSeparationAndForce(double &,double &,double,double,double,double) const;
		virtual double GetSslideAcDt(double,double,double,double,bool &,double) const;
	
		// accessors
		virtual const char *MaterialType(void) const;
		void SetFrictionCoeff(double);
		virtual bool IgnoreContact(void) const;
		virtual bool ContactLawNeedsContactArea(void) const;
		virtual bool IsFrictionless(void) const;
		virtual bool IsStick(void) const;
	
	protected:
		double frictionCoeff;
		double frictionCoeffStatic;
		int frictionStyle;
		double displacementOnly;
		double Dc;
	
};

#endif
