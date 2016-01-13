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
		LinearInterface();
		LinearInterface(char *matName);

		// initialize
		virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintContactLaw(void) const;
	
		// methods
		virtual bool GetInterfaceForcesForNode(Vector *,Vector *,double *,double,Vector *,double,bool,bool,double,
										   Vector *,double,double) const;
		virtual bool GetCrackInterfaceForce(Vector *,Vector *,Vector *,double,double,Vector *,double *) const;

		// accessors
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag() const;
		void SetParameters(double,double,double);
		virtual bool IgnoreContact(void) const;
		virtual bool IsImperfectInterface(void) const;
		virtual bool IsPerfectTangentialInterface(void) const;
		virtual bool IsPerfectNormalInterface(bool) const;
		virtual bool IsPerfectNormalInterface(void) const;
	
	protected:
		double Dnt,Dnc,Dt;
		bool hasSetDnc;
	
};

#endif
