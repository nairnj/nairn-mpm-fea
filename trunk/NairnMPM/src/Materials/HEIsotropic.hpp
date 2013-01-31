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

#include "Materials/HyperElastic.hpp"

class HardeningLawBase;

//enum {G1_PROP=0,G2_PROP=0,KBULK_PROP,CTEE_PROP,HEISOTROPIC_PROPS};

class HEIsotropic : public HyperElastic
{
    public:
        double G1,G2;
		// double aI,betaI		// isotropic expanion defined in super classes
		// JAN: never used as class variable
        //double magnitude_strial;
		// JAN: deleted Ep, yield, gyld. First two in hardening law; gyld not used

        // constructors and destructors
		HEIsotropic();
		HEIsotropic(char *matName);
		
		// initialize
        virtual char *InputMat(char *,int &);
        virtual const char *VerifyProperties(int);
        virtual void ValidateForUse(int);
		virtual void PrintMechanicalProperties(void);
		virtual void InitialLoadMechProps(int,int);
        virtual char *InitHistoryData(void);
		
		// step methods
		virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,int);
		virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int);
    
        Tensor GetTrialStressTensor2D(Tensor *,double);
        Tensor GetTrialStressTensor3D(Tensor *,double);
        Tensor GetNormalTensor2D(Tensor *,double);
        Tensor GetNormalTensor3D(Tensor *,double);
		virtual double GetDilationalTerms(MPMBase *,double,int,double &);
    
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
		virtual double WaveSpeed(bool,MPMBase *);
        virtual double GetHistory(int,char *);
        virtual double GetMagnitudeS(Tensor *st,int);	
    
    protected:
		// JAN: deleted yield properties, dalpha, and alpint
        // unique properties
        double Gred,Kred;
        double G1sp,G2sp;

		// JAN: for future new plastic laws, currently hard-coded to linear law
		HardeningLawBase *plasticLaw;
	
		// JAN: J history might move depening on hardening law
		int J_history;
};

#endif
