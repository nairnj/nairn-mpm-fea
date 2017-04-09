/********************************************************************************
	NoninearInterface.hpp
	nairn-mpm-fea
 
	Created by John Nairn, 3/30/2017.
	Copyright (c) 2017 John A. Nairn, All rights reserved.
 
	Dependencies
		LinearInterface.hpp,ContactLaw.hpp,MaterialBase.hpp
********************************************************************************/

#ifndef NONLINEARINTERFACELAW

#define NONLINEARINTERFACELAW 65

#include "Materials/LinearInterface.hpp"

#define FORCE_STICK -1
#define STABLE 0
#define FORCE_DEBOND 1

class NonlinearInterface : public LinearInterface
{
	public:
		
		// constructors and destructors
		NonlinearInterface();
		NonlinearInterface(char *matName);
		
		// initialize
		virtual char *InputMaterialProperty(char *,int &,double &);
		virtual void PrintContactLaw(void) const;
		
		// methods
		virtual void GetInterfaceForces(Vector *,Vector *,double *,double,Vector *,double,double,Vector *,double,double,double) const;
	
		// interface law to implement
		virtual int CheckDtStability(double,double) const;
		virtual double GetFt(double,double) const;
		virtual double GetFtPrime(double,double) const;
		virtual double GetFtEnergy(double,double) const;
	
		virtual int CheckDnStability(double,double,double) const;
		virtual double GetFn(double,double) const;
		virtual double GetFnPrime(double,double) const;
		virtual double GetFnEnergy(double,double) const;
	
		// accessors
		virtual const char *MaterialType(void) const;
	
	protected:
		int order;
	
};

#endif
