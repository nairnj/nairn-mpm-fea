/********************************************************************************
    HardeningLawBase.hpp
    NairnMPM

    Created by John Nairn, 1/17/2103
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Dependencies
        MaterialBase.hpp
********************************************************************************/

#ifndef _HARDENINGLAWBASE_

#define _HARDENINGLAWBASE_

#include "Materials/MaterialBase.hpp"

class MPMBase;

// plastic law properties
typedef struct {
	double alpint;
	double dalpha;
} HardeningAlpha;

class HardeningLawBase
{
    public:
        HardeningLawBase();
        HardeningLawBase(MaterialBase *);
        virtual ~HardeningLawBase();
    
        // initialize
        virtual char *InputMat(char *,int &);
		virtual const char *VerifyAndLoadProperties(int);
        virtual void PrintYieldProperties(void) const = 0;
		virtual int HistoryDoublesNeeded(void) const;
    
        // load properties
		virtual void *GetCopyOfHardeningProps(MPMBase *,int);
		virtual void DeleteCopyOfHardeningProps(void *,int) const;
        virtual double GetShearRatio(MPMBase *,double,double,void *) const;
	
		// hardening law core methods
        virtual double GetYield(MPMBase *,int,double,HardeningAlpha *,void *) const = 0;
        virtual double GetYieldIncrement(MPMBase *,int,double,HardeningAlpha *,void *) const;
        virtual double GetKPrime(MPMBase *,int,double,HardeningAlpha *,void *) const = 0;
        virtual double GetK2Prime(MPMBase *,double,double,HardeningAlpha *,void *) const = 0;
    
        // return mapping methods
        virtual double SolveForLambda(MPMBase *,int,double,Tensor *,double,double,double,double,HardeningAlpha *,void *) const;
        virtual double SolveForLambdaBracketed(MPMBase *,int,double,Tensor *,double,double,double,double,HardeningAlpha *,void *) const;
        virtual bool LambdaConverged(int,double,double) const;
	
        // Default internal variable as cumlative plastic strain
        virtual void UpdateTrialAlpha(MPMBase *,int,HardeningAlpha *) const;
        virtual void UpdateTrialAlpha(MPMBase *,int,double,double,HardeningAlpha *) const;
        virtual void UpdatePlasticInternal(MPMBase *,int,HardeningAlpha *) const;
		virtual void ElasticUpdateFinished(MPMBase *,int,double) const;
	
        // accessors
		virtual double GetHistory(int,char *) const;
        virtual const char *GetHardeningLawName(void) const = 0;
    
    protected:
        double yield,yldred;
        MaterialBase *parent;
	
		virtual void BracketSolution(MPMBase *,int,double,Tensor *,double,double,double,double,double *,double *,HardeningAlpha *a,void *) const;
};

#endif
