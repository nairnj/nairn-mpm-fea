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

class HardeningLawBase
{
    public:
        HardeningLawBase();
        HardeningLawBase(MaterialBase *);
        virtual ~HardeningLawBase();
    
        // initialize
        virtual char *InputMat(char *,int &);
        virtual void PrintYieldProperties(void) = 0;
        virtual void InitialLoadMechProps(int,int);
        virtual const char *VerifyProperties(int);
		virtual int HistoryDoublesNeeded(void);
    
        // hardening law core methods
        virtual double GetShearRatio(MPMBase *,double,double);
        virtual void LoadHardeningLawProps(MPMBase *,int);
        virtual double GetYield(MPMBase *,int,double) = 0;
        virtual double GetYieldIncrement(MPMBase *,int,double);
        virtual double GetKPrime(MPMBase *,int,double) = 0;
        virtual double GetK2Prime(MPMBase *,double,double) = 0;
    
        // return mapping methods
        virtual double SolveForLambda(MPMBase *,int,double,Tensor *,double,double,double,double);
        virtual double SolveForLambdaBracketed(MPMBase *,int,double,Tensor *,double,double,double,double);
        virtual void BracketSolution(MPMBase *,int,double,Tensor *,double,double,double,double,double *,double *);
        virtual bool LambdaConverged(int,double,double);
	
        // Default internal variable as cumlative plastic strain
        virtual void UpdateTrialAlpha(MPMBase *,int);
        virtual void UpdateTrialAlpha(MPMBase *,int,double,double);
        virtual void UpdatePlasticInternal(MPMBase *,int);
		virtual void ElasticUpdateFinished(MPMBase *,int,double);
	
        // accessors
		virtual double GetHistory(int,char *);
        virtual const char *GetHardeningLawName(void) = 0;
    
    protected:
        double yield,yldred;
        MaterialBase *parent;
        double alpint,dalpha;
};

#endif
