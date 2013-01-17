/********************************************************************************
 MGJohnsonCook.hpp
 NairnMPM
 
 Created by John Nairn, 1/17/13.
 Copyright (c) 2012 John A. Nairn, All rights reserved.
 
 Dependencies
    MGSCGLMaterial.hpp (Isoplasticity.hpp IsotropicMat.hpp,
                Elastic.hpp, MaterialBase.hp)
 ********************************************************************************/

#ifndef MGJOHNSONCOOK

#define MGJOHNSONCOOK 25

#include "Materials/MGSCGLMaterial.hpp"

class MGJohnsonCook : public MGSCGLMaterial
{
    public:
        
        // constructors and destructors
        MGJohnsonCook();
        MGJohnsonCook(char *);
        
        // initialize
        virtual char *InputMat(char *,int &);
        virtual void PrintYieldProperties(void);
        virtual void InitialLoadMechProps(int,int);
        
        // methods
        virtual void LoadMechanicalProps(MPMBase *,int);
        
        // override plastic potential functions
        virtual double GetYield(MPMBase *,int,double);
        virtual double GetKPrime(MPMBase *,int,double);
        virtual double GetK2Prime(MPMBase *,double,double);
        
        // accessors
        virtual const char *MaterialType(void);
        virtual int MaterialTag();
        
    protected:
        double Bjc,Cjc,njc,ep0jc,Tmjc,mjc;
        double Bred,TjcTerm,edotMin,eminTerm;
        
};

#endif
