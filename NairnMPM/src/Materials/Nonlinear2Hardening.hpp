/********************************************************************************
	Nonlinear2Hardening.hpp
	NairnMPM

	Created by John Nairn, 2/8/2103
	Copyright (c) 2013 John A. Nairn, All rights reserved.

	Dependencies
		NonlinearHardening.hpp, HardeningLawBase.hpp, MaterialBase.hpp
********************************************************************************/

#ifndef _NONLINEAR2HARDENING_

#define _NONLINEAR2HARDENING_

#include "Materials/NonlinearHardening.hpp"

class Nonlinear2Hardening : public NonlinearHardening
{
	public:
		Nonlinear2Hardening();
		Nonlinear2Hardening(MaterialBase *);
    
		// hardening law core methods
		virtual double GetYield(MPMBase *,int,double);
		virtual double GetKPrime(MPMBase *,int,double);
		virtual double GetK2Prime(MPMBase *,double,double);
    
		// accessors
		virtual const char *GetHardeningLawName(void);
};

#endif
