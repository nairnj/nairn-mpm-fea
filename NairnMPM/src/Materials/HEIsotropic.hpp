/********************************************************************************
    HEIsotropic.hpp
    nairn-mpm-fea
    
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
        // constructors and destructors
		HEIsotropic(char *,int);
		
		// initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
        virtual bool AcceptHardeningLaw(HardeningLawBase *,int );
		virtual HardeningLawBase *GetPlasticLaw(void) const;
		virtual const char *VerifyAndLoadProperties(int);
	
		// history data
		virtual int SizeOfHistoryData(void) const;
   		virtual int NumberOfHistoryDoubles(void) const;
		virtual char *InitHistoryData(char *,MPMBase *);
	
		// const methods
        virtual void ValidateForUse(int) const;
		virtual void PrintMechanicalProperties(void) const;
		
		// step methods
        virtual int SizeOfMechanicalProperties(int &) const;
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *,int) const;
	
		// constitutive law
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int) const;
        virtual void UpdatePressure(MPMBase *,double,double,int,double,double,HEPlasticProperties *,ResidualStrains *,
									double,int,double &,double &) const;
        Tensor GetTrialDevStressTensor(Tensor *,double,int,double) const;
        virtual double GetMagnitudeS(Tensor *st,int) const;
        Tensor GetNormalTensor(Tensor *,double,int) const;
    
		// accessors
        virtual Vector ConvertJToK(Vector,Vector,Vector,int);
        virtual Tensor GetStress(Tensor *,double,MPMBase *) const;
		virtual void SetStress(Tensor *,MPMBase *) const;
		virtual void IncrementThicknessStress(double,MPMBase *) const;
		virtual const char *MaterialType(void) const;
		virtual double WaveSpeed(bool,MPMBase *) const;
        virtual double CurrentWaveSpeed(bool,MPMBase *,int) const;
        virtual bool SupportsArtificialViscosity(void) const;
		virtual int AltStrainContains(void) const;
    
    protected:
		double G1,gammaI;
		// double aI,betaI		// isotropic expansion defined in super classes
		double G1sp;
		HardeningLawBase *plasticLaw;
};

#endif

