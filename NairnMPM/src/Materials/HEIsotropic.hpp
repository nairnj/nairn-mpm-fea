/********************************************************************************
    HEIsotropic.hpp
    NairnMPM
    
    Created by John Nairn, Sept 27, 2011.
    Copyright (c) 2011 John A. Nairn, All rights reserved.

	Dependencies
		HyperElastic.hpp (MaterialBase.hpp)
********************************************************************************/

#ifndef HEISOTROPIC

#define HEISOTROPIC 24

#define SQRT_TWOTHIRDS 0.8164965809277260
#define TWOTHIRDS 0.6666666666666667
#define ONETHIRD 0.3333333333333333

#define J_HISTORY 0
#define ALPHA_HISTORY 1


#include "Materials/HyperElastic.hpp"

//enum {G1_PROP=0,G2_PROP=0,KBULK_PROP,CTEE_PROP,HEISOTROPIC_PROPS};

class HEIsotropic : public HyperElastic
{
    public:
        double G1,G2;
        double magnitude_strial;
        // double aI,betaI		// isotropic expanion defined in super classes
   
        // Plastic modulus (Ep: slope of unidirectional stress - plastic strain curve) or tangential modulus ET
        // can enter Ep OR ET for linear hardening or enter beta AND npow for nonlinear hardening

        double Ep;
        double yield, gyld;

        // constructors and destructors
		HEIsotropic();
		HEIsotropic(char *matName);
		
		// initialize
        virtual char *InputMat(char *,int &);
        virtual const char *VerifyProperties(int);
        virtual void ValidateForUse(int);
		virtual void PrintMechanicalProperties(void);
		virtual void InitialLoadMechProps(int,int);
        virtual char *MaterialData(void);
		
		// step methods
		virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,int);
		virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int);
    
        //virtual double GetYield(MPMBase *,int,double,double);
        virtual double GetYield(MPMBase *,int,double);
        //virtual double SolveForLambda(MPMBase *,int,double,Tensor *,double);
        
        Tensor GetTrialStressTensor2D(Tensor *,double);
        Tensor GetTrialStressTensor3D(Tensor *,double);
        Tensor GetNormalTensor2D(Tensor *,double);
        Tensor GetNormalTensor3D(Tensor *,double);
        //Tensor GetMagnitudeS(double,int);
    
    
        // Default internal variable as cumlative plastic strain
        virtual void UpdateTrialAlpha(MPMBase *,int);
        virtual void UpdateTrialAlpha(MPMBase *,int,double,double);
        virtual void UpdatePlasticInternal(MPMBase *,int);
				
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
		virtual double WaveSpeed(bool,MPMBase *);
        virtual double GetHistory(int,char *);
        virtual double GetMagnitudeS(Tensor *st,int);	
    
    protected:
        // unique properties
        double yldred,Gred,Kred;
        double Epred;
        double alpint,dalpha;		// internal alpha variable and current increment
        double G1sp, G2sp;

};

#endif

