/********************************************************************************
    SCGLHardening.hpp
    NairnMPM

    Created by John Nairn, 1/18/2103
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Dependencies
        HardeningLawBase.hpp, MaterialBase.hpp
********************************************************************************/

#ifndef _SCGLHARDENING_

#define _SCGLHARDENING_

#include "Materials/HardeningLawBase.hpp"

class SCGLHardening : public HardeningLawBase
{
    public:
        // contructors
        SCGLHardening();
        SCGLHardening(MaterialBase *);
        
        // initialize
        virtual char *InputMat(char *,int &);
        virtual const char *VerifyProperties(int);
        virtual void PrintYieldProperties(void);
        virtual void InitialLoadMechProps(int,int);
        
        // hardening law core methods
        virtual double GetShearRatio(MPMBase *,double,double);
        virtual double GetYield(MPMBase *,int,double);
        virtual double GetKPrime(MPMBase *,int,double);
        virtual double GetK2Prime(MPMBase *,double,double);
        virtual double GetYieldIncrement(MPMBase *,int,double);
        
        // accessors
        virtual const char *GetHardeningLawName(void);
    
    protected:
        double GPp,GTp,beta,nhard,yieldMax;
        double GPpred,yldMaxred,Gratio;
};

#endif
