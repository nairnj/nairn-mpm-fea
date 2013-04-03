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

// plastic law properties
typedef struct {
	double Gred;
	double Kred;
	void *hardProps;
} HEPlasticProperties;

class HardeningLawBase;

//enum {G1_PROP=0,G2_PROP=0,KBULK_PROP,CTEE_PROP,HEISOTROPIC_PROPS};

class HEIsotropic : public HyperElastic
{
    public:
        double G1;
		// double aI,betaI		// isotropic expansion defined in super classes

        // constructors and destructors
		HEIsotropic();
		HEIsotropic(char *matName);
		
		// initialize
        virtual char *InputMat(char *,int &);
        virtual void SetHardeningLaw(char *);
        virtual const char *VerifyAndLoadProperties(int);
		virtual char *InitHistoryData(void);
	
		// const methods
        virtual void ValidateForUse(int) const;
		virtual void PrintMechanicalProperties(void) const;
		
		// step methods
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int);
		virtual void DeleteCopyOfMechanicalProps(void *,int) const;
	
		// constitutive law
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *);
        virtual void UpdatePressure(MPMBase *,double,double,int,double,double,HEPlasticProperties *,ResidualStrains *);
        Tensor GetTrialDevStressTensor(Tensor *,double,int,double) const;
        virtual double GetMagnitudeS(Tensor *st,int) const;
        Tensor GetNormalTensor(Tensor *,double,int) const;
    
		// accessors
        virtual Tensor GetStress(Tensor *,double) const;
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag() const;
		virtual double WaveSpeed(bool,MPMBase *) const;
        virtual double GetHistory(int,char *) const;
        virtual bool SupportsArtificialViscosity(void) const;
    
    protected:
        double G1sp;
		HardeningLawBase *plasticLaw;
		int J_history;						// History variable may move depending on plastic law
};

#endif

