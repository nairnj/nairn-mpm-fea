/********************************************************************************
    LinearHardening.hpp
    NairnMPM

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
        virtual void PrintYieldProperties(void);
        virtual void InitialLoadMechProps(int,int);
    
        // hardening law core methods
        virtual double GetYield(MPMBase *,int,double);
        virtual double GetKPrime(MPMBase *,int,double);
        virtual double GetK2Prime(MPMBase *,double,double);
    
        // return mapping can have fast option
        virtual double SolveForLambdaBracketed(MPMBase *,int,double,Tensor *,double,double,double);
    
		// hyperelastic return mapping methods
		virtual double HESolveForLambdaBracketed(MPMBase *,int,double,double,double);
    
        // accessors
        virtual const char *GetHardeningLawName(void);
    
    protected:
        double beta,Epred;
};

#endif
