/********************************************************************************
    JohnsonCook.hpp
    NairnMPM
    
    Created by John Nairn, August 12, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		HardeningLawBase.hpp, MaterialBase.hpp
********************************************************************************/

#ifndef _JOHNSONCOOKHARDENING_

#define _JOHNSONCOOKHARDENING_

#include "Materials/HardeningLawBase.hpp"

class JohnsonCook : public HardeningLawBase
{
    public:
        // contructors
        JohnsonCook();
        JohnsonCook(MaterialBase *);
        
        // initialize
        virtual char *InputMat(char *,int &);
        virtual void PrintYieldProperties(void);
        virtual void InitialLoadMechProps(int,int);
        
        // hardening law core methods
        virtual void LoadHardeningLawProps(MPMBase *,int);
        virtual double GetYield(MPMBase *,int,double);
        virtual double GetKPrime(MPMBase *,int,double);
        virtual double GetK2Prime(MPMBase *,double,double);
        virtual double GetYieldIncrement(MPMBase *,int,double);
    
        // accessors
        virtual const char *GetHardeningLawName(void);
    
    protected:
		double Bjc,Cjc,njc,ep0jc,Tmjc,mjc;
		double Bred,TjcTerm,edotMin,eminTerm;

};

#endif
