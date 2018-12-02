/********************************************************************************
	LinearInterface.hpp
	nairn-mpm-fea

	Created by John Nairn, Oct 24, 2015.
	Copyright (c) 2015 John A. Nairn, All rights reserved.

	Dependencies
		ContactLaw.hpp,MaterialBase.hpp
********************************************************************************/

#ifndef LINEARINTERFACELAW

#define LINEARINTERFACELAW 62

#include "Materials/ContactLaw.hpp"

class LinearInterface : public ContactLaw
{
	public:

		// constructors and destructors
		LinearInterface(char *,int);

		// initialize
		virtual char *InputContactProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintContactLaw(void) const;
	
		// methods
		virtual void GetInterfaceForces(Vector *,Vector *,double *,double,Vector *,double,double,Vector *,double,double,double) const;
	
		// accessors
		virtual const char *MaterialType(void) const;
		void SetParameters(double,double,double);
		virtual bool IgnoreContact(void) const;
		virtual bool IsImperfectInterface(void) const;
	
	protected:
		double Dnt,Dnc,Dt;
		bool hasSetDnc;					// false means linear in normal direction
	
};

#endif
