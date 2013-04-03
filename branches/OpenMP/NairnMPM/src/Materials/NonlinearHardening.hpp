/********************************************************************************
    NonlinearHardening.hpp
    NairnMPM

    Created by John Nairn, 1/17/2103
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Dependencies
        HardeningLawBase.hpp, MaterialBase.hpp
********************************************************************************/

#ifndef _NONLINEARHARDENING_

#define _NONLINEARHARDENING_

#include "Materials/HardeningLawBase.hpp"

class NonlinearHardening : public HardeningLawBase
{
    public:
        NonlinearHardening();
        NonlinearHardening(MaterialBase *);
    
        virtual char *InputMat(char *,int &);
        virtual void PrintYieldProperties(void);
    
        // hardening law core methods
        virtual double GetYield(MPMBase *,int,double);
        virtual double GetKPrime(MPMBase *,int,double);
        virtual double GetK2Prime(MPMBase *,double,double);
    
        // accessors
        virtual const char *GetHardeningLawName(void);
    
    protected:
        double beta,npow;
};

#endif
