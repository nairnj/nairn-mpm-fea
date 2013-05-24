/********************************************************************************
	Nonlinear2Hardening.hpp
	nairn-mpm-fea

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
		virtual double GetYield(MPMBase *,int,double,HardeningAlpha *,void *) const;
		virtual double GetKPrime(MPMBase *,int,double,HardeningAlpha *,void *) const;
		virtual double GetK2Prime(MPMBase *,double,double,HardeningAlpha *,void *) const;
    
		// accessors
		virtual const char *GetHardeningLawName(void) const;
};

#endif
