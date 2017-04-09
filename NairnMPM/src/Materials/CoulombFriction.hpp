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
		CoulombFriction();
		CoulombFriction(char *matName);
	
		// initialize
		virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintContactLaw(void) const;
	
		// methods
		virtual bool GetFrictionalDeltaMomentum(Vector *,Vector *,double,double *,double,bool,double,bool,double,Vector *) const;
		virtual double GetSslideAcDt(double,double,double,double,bool &,double) const;
	
		// accessors
		virtual const char *MaterialType(void) const;
		void SetFrictionCoeff(double);
		virtual bool IgnoreContact(void) const;
	
	protected:
		double frictionCoeff;
		double frictionCoeffStatic;
		int frictionStyle;
	
};

#endif
