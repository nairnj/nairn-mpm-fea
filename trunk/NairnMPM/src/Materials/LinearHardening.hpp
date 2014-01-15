/********************************************************************************
    LinearHardening.hpp
    nairn-mpm-fea

    Created by John Nairn, 1/17/2103
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Dependencies
        HardeningLawBase.hpp, MaterialBase.hpp
********************************************************************************/

#ifndef _LINEARHARDENING_

#define _LINEARHARDENING_

#include "Materials/HardeningLawBase.hpp"

class LinearHardening : public HardeningLawBase
{
    public:
        double Ep;
    
        // contructors
        LinearHardening();
        LinearHardening(MaterialBase *);
    
        // initialize
        virtual char *InputMat(char *,int &);
        virtual void PrintYieldProperties(void) const;
		virtual const char *VerifyAndLoadProperties(int);
    
        // hardening law core methods
        virtual double GetYield(MPMBase *,int,double,HardeningAlpha *,void *) const;
        virtual double GetYieldIncrement(MPMBase *,int,double,HardeningAlpha *,void *) const;
        virtual double GetKPrime(MPMBase *,int,double,HardeningAlpha *,void *) const;
        virtual double GetK2Prime(MPMBase *,double,double,HardeningAlpha *,void *) const;
    
        // return mapping can have fast option
        virtual double SolveForLambdaBracketed(MPMBase *,int,double,Tensor *,double,double,double,double,HardeningAlpha *,void *) const;
    
        // accessors
        virtual const char *GetHardeningLawName(void) const;
    
    protected:
        double beta,Epred;
};

#endif
