/********************************************************************************
    VonMisesHardening.hpp
    NairnMPM
    
    Created by Yajun Guo in Jan 2005.
    Copyright (c) 2005 John A. Nairn, All rights reserved.

	Dependencies
		Isoplasticity.hpp (IsotropicMat.hpp, Elastic.hpp, MaterialBase.hp)
********************************************************************************/

#ifndef VONMISESHARDENING

#define VONMISESHARDENING 9

#include "Materials/IsoPlasticity.hpp"

class VonMisesHardening : public IsoPlasticity
{
    public:
        // Plastic modulus (Ep: slope of unidirectional stress - plastic strain curve) or tangential modulus ET
		// yield stress in IsoPlasticity
        double Ep,ET;

        // constructors and destructors
        VonMisesHardening();
        VonMisesHardening(char *);
        
        // initialize
        virtual char *InputMat(char *,int &);
		virtual void InitialLoadMechProps(int,int);
		virtual void PrintMechanicalProperties(void);
		virtual void PrintYieldProperties(void);
				
		// override plastic potential functions
		virtual double GetYield(MPMBase *,int,double);
 		virtual double GetKPrime(MPMBase *,int,double);
		virtual double GetK2Prime(MPMBase *,double,double);
		virtual double SolveForLambda(MPMBase *,int,double,Tensor *,double);
       
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
		
    protected:
        // specific values of yield strength & plastic modulus 
        double Epred;

};

#endif
