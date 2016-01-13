/********************************************************************************
    JohnsonCook.hpp
    nairn-mpm-fea
    
    Created by John Nairn, August 12, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		HardeningLawBase.hpp, MaterialBase.hpp
********************************************************************************/

#ifndef _JOHNSONCOOKHARDENING_

#define _JOHNSONCOOKHARDENING_

#include "Materials/HardeningLawBase.hpp"

// plastic law properties
typedef struct {
	double hmlgTemp;
	double TjcTerm;
} JCProperties;

class JohnsonCook : public HardeningLawBase
{
    public:
        // contructors
        JohnsonCook();
        JohnsonCook(MaterialBase *);
        
        // initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
        virtual void PrintYieldProperties(void) const;
		virtual const char *VerifyAndLoadProperties(int);
	
		// copy of properties
        virtual int SizeOfHardeningProps(void) const;
		virtual void *GetCopyOfHardeningProps(MPMBase *,int,void *,int);
		virtual void DeleteCopyOfHardeningProps(void *,int) const;
    
        // hardening law core methods
        virtual double GetYield(MPMBase *,int,double,HardeningAlpha *,void *) const;
        virtual double GetKPrime(MPMBase *,int,double,HardeningAlpha *,void *) const;
        virtual double GetK2Prime(MPMBase *,double,double,HardeningAlpha *,void *) const;
        virtual double GetYieldIncrement(MPMBase *,int,double,HardeningAlpha *,void *) const;
	
		// return mapping
        virtual double SolveForLambdaBracketed(MPMBase *,int,double,Tensor *,double,double,
											   double,double,HardeningAlpha *a,void *,int) const;
    
        // accessors
        virtual const char *GetHardeningLawName(void) const;
    
    protected:
		double Bjc,Cjc,njc,ep0jc,Tmjc,mjc;
		double Bred,edotMin,eminTerm;

};

#endif
