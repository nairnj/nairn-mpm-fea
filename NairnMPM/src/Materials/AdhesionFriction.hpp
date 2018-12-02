/********************************************************************************
	AdhesionFriction.hpp
	nairn-mpm-fea

	Created by John Nairn, Nov 16, 2015.
	Copyright (c) 2015 John A. Nairn, All rights reserved.

	Dependencies
		CoulombFriction.hpp,ContactLaw.hpp,MaterialBase.hpp
********************************************************************************/

#ifndef ADHESIONFRICTIONLAW

#define ADHESIONFRICTIONLAW 63

#include "Materials/CoulombFriction.hpp"

class AdhesionFriction : public CoulombFriction
{
	public:
		
		// constructors and destructors
		AdhesionFriction(char *,int);
		
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
		double Sa;
		double Na;
		double kmu;
		double vhalf;
		bool smoothStaticToDynamic;
	
};

#endif
