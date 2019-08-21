/********************************************************************************
    HardeningLawBase.hpp
    nairn-mpm-fea

    Created by John Nairn, 1/17/2103
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Dependencies
        MaterialBase.hpp
********************************************************************************/

#ifndef _HARDENINGLAWBASE_

#define _HARDENINGLAWBASE_

#include "Materials/MaterialBase.hpp"

class MPMBase;

class HardeningLawBase
{
    public:
    
		HardeningLawBase();
        HardeningLawBase(MaterialBase *);
        virtual ~HardeningLawBase();
    
        // initialize
        virtual char *InputHardeningProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
        virtual void PrintYieldProperties(void) const = 0;
	
		// history data
		virtual int HistoryDoublesNeeded(void) const;
        virtual void InitPlasticHistoryData(double *) const;
		virtual double GetHistory(int,char *) const;
    
        // load properties
        virtual int SizeOfHardeningProps(void) const;
		virtual void *GetCopyOfHardeningProps(MPMBase *,int,void *,int);
		virtual void DeleteCopyOfHardeningProps(void *,int) const;
        virtual double GetShearRatio(MPMBase *,double,double,void *,int) const;
	
		// hardening law core methods
        virtual double GetYield(MPMBase *,int,double,HardeningAlpha *,void *) const = 0;
        virtual double GetYieldIncrement(MPMBase *,int,double,HardeningAlpha *,void *) const;
        virtual double GetKPrime(MPMBase *,int,double,HardeningAlpha *,void *) const = 0;
        virtual double GetK2Prime(MPMBase *,double,double,HardeningAlpha *,void *) const = 0;
    
        // return mapping methods
        virtual double SolveForLambda(MPMBase *,int,double,Tensor *,double,
									  double,double,double,HardeningAlpha *,void *,int) const;
        virtual double SolveForLambdaBracketed(MPMBase *,int,double,Tensor *,double,
											   double,double,double,HardeningAlpha *,void *,int) const;
        virtual bool LambdaConverged(int,double,double) const;
	
        // Default internal variable as cumlative plastic strain
        virtual void UpdateTrialAlpha(MPMBase *,int,HardeningAlpha *,int) const;
        virtual void UpdateTrialAlpha(MPMBase *,int,double,double,HardeningAlpha *,int) const;
        virtual void UpdatePlasticInternal(MPMBase *,int,HardeningAlpha *,int) const;
		virtual void ElasticUpdateFinished(MPMBase *,int,double,int) const;
	
        // accessors
        virtual const char *GetHardeningLawName(void) const = 0;
		virtual int GetHardeningID(void) const;
    
    protected:
        double yield,yldred;
		double yieldMin,yldredMin;
        MaterialBase *parent;
		int lawID;
	
		virtual void BracketSolution(MPMBase *,int,double,Tensor *,double,double,double,
									 double,double *,double *,HardeningAlpha *a,void *,int) const;
};

#endif
