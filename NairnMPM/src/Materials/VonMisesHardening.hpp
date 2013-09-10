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

#define YT_HISTORY 1
#define EPDOT_HISTORY 2

#include "Materials/IsoPlasticity.hpp"

class VonMisesHardening : public IsoPlasticity
{
    public:
        // Plastic modulus (Ep: slope of unidirectional stress - plastic strain curve) or tangential modulus ET
		// can enter Ep OR ET for linear hardening or enter beta AND npow for nonlinear hardening
		// yield stress in IsoPlasticity
        double Ep,ET;
		double beta,npow;

        // constructors and destructors
        VonMisesHardening();
        VonMisesHardening(char *);
        
        // initialize
        virtual char *InputMat(char *,int &);
		virtual void InitialLoadMechProps(int,int);
		virtual void PrintMechanicalProperties(void);
		virtual void PrintYieldProperties(void);
		virtual char *MaterialData(void);
				
		// override plastic potential functions
		virtual double GetYield(MPMBase *,int,double);
 		virtual double GetKPrime(MPMBase *,int,double);
		virtual double GetK2Prime(MPMBase *,double,double);
		virtual double SolveForLambdaBracketed(MPMBase *,int,double,Tensor *,double);
		virtual void UpdatePlasticInternal(MPMBase *,int);
       
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
		virtual double GetHistory(int,char *);
		
    protected:
        // specific values of yield strength & plastic modulus 
        double Epred;
		bool linearHardening;

};

#endif
