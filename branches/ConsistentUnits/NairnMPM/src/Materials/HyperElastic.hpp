/********************************************************************************
    HyperElastic.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 24 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
 
    Base class for all Hyperelastic materials. It should have methods common
    to all.

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef _HYPERELASTIC_

#define _HYPERELASTIC_

#include "Materials/MaterialBase.hpp"

enum { HALF_J_SQUARED_MINUS_1_MINUS_LN_J=0,J_MINUS_1_SQUARED,LN_J_SQUARED };

class HyperElastic : public MaterialBase
{
    public:
	
        // constructors and destructors
        HyperElastic();
        HyperElastic(char *);
    
        // initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
        virtual void SetInitialParticleState(MPMBase *,int) const;
    
		// Methods (make virtual if any subclass needs them)
        double IncrementDeformation(MPMBase *,Matrix3,Tensor *,int) const;
		double GetResidualStretch(MPMBase *,double &,ResidualStrains *) const;
        virtual double GetCpMinusCv(MPMBase *) const;
    
        // Accessors
        virtual double GetVolumetricTerms(double,double) const;
        virtual void GetNewtonPressureTerms(double,double,double &,double &) const;
        virtual double GetCurrentRelativeVolume(MPMBase *) const;
    
    protected:
		double Kbulk;               // bulk modulus
		double aI;                  // thermal expansion isotropic
		// double betaI;			// moisture expansion isotopic (in base material)
	
        int UofJOption;             // pick U(J) function
        double Ksp;                 // specific bulk modulus
        double Ka2sp;               // For Cp-Cv
};

#endif

