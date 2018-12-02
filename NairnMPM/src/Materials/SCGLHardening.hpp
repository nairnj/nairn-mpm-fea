/********************************************************************************
    SCGLHardening.hpp
    nairn-mpm-fea

    Created by John Nairn, 1/18/2103
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Dependencies
        HardeningLawBase.hpp, MaterialBase.hpp
********************************************************************************/

#ifndef _SCGLHARDENING_

#define _SCGLHARDENING_

#include "Materials/HardeningLawBase.hpp"

// plastic law properties
typedef struct {
	double Gratio;
} SCGLProperties;

class SCGLHardening : public HardeningLawBase
{
    public:
        // contructors
        SCGLHardening();
        SCGLHardening(MaterialBase *);
        
        // initialize
        virtual char *InputHardeningProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
        virtual void PrintYieldProperties(void) const;
    
		// copy of properties
        virtual int SizeOfHardeningProps(void) const;
		virtual void *GetCopyOfHardeningProps(MPMBase *,int,void *,int);
		virtual void DeleteCopyOfHardeningProps(void *,int) const;
		virtual double GetShearRatio(MPMBase *,double,double,void *,int) const;
	
        // hardening law core methods
        virtual double GetYield(MPMBase *,int,double,HardeningAlpha *a,void *) const;
        virtual double GetKPrime(MPMBase *,int,double,HardeningAlpha *a,void *) const;
        virtual double GetK2Prime(MPMBase *,double,double,HardeningAlpha *a,void *) const;
        virtual double GetYieldIncrement(MPMBase *,int,double,HardeningAlpha *,void *) const;
        
        // accessors
        virtual const char *GetHardeningLawName(void) const;
    
    protected:
		// properties independent of particle state
        double GPp,GTp,beta,nhard,yieldMax;
		double GPpred,yldMaxred;
};

#endif
