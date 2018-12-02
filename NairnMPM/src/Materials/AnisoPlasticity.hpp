/********************************************************************************
    AnisoPlasticity.hpp
    nairn-mpm-fea
    
    Created by John Nairn, June 12, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		Orthotropic.hpp (TranIsotropic.hpp, Elastic.hpp, MaterialBase.hpp)
********************************************************************************/

#ifndef _ANISOPLASTICITY_

#define _ANISOPLASTICITY_

#include "Materials/Orthotropic.hpp"

// plastic law properties
typedef struct {
	double aint;
	double sAsmag;
	Tensor dfds;
	double dfCdf;
	Tensor Cdf;
	ElasticProperties *ep;
} AnisoPlasticProperties;

class AnisoPlasticity : public Orthotropic
{
    public:
		static int warnNonconvergence;
    
        // constructors and destructors
		AnisoPlasticity(char *,int);
		
        // initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void ValidateForUse(int) const;
		virtual void PrintMechanicalProperties(void) const;
		virtual void PrintYieldProperties(void) const;
			
		// methods
        virtual int SizeOfMechanicalProperties(int &) const;
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *,int) const;
		virtual int AltStrainContains(void) const;
		virtual ElasticProperties *GetElasticPropertiesPointer(void *) const;
	
		// Large Rotation Methods
		virtual void LRElasticConstitutiveLaw(MPMBase *,Matrix3 &,Matrix3 &,Matrix3 &,Matrix3 &,Matrix3 &,int,void *,ResidualStrains *) const;
		virtual Tensor SolveForPlasticIncrement(MPMBase *,int,double,Tensor &,AnisoPlasticProperties *) const;
		virtual double GetMagnitudeHill(Tensor &,int) const;
		virtual void GetDfDsigma(Tensor &,int,AnisoPlasticProperties *) const;
		virtual void GetDfCdf(Tensor &,int,AnisoPlasticProperties *) const;
	
		// hardening term methods (move to hardening law class when want more hardening options)
		virtual double GetYield(AnisoPlasticProperties *p) const = 0;
		virtual double GetGPrime(AnisoPlasticProperties *) const = 0;
		
   protected:
		double syxx,syyy,syzz,tyyz,tyxz,tyxy;
		double syxxred2,syyyred2,syzzred2,tyyzred2,tyxzred2,tyxyred2;		// equal to 1/yield^2 and reduced
		double fTerm,gTerm,hTerm,sigmaYref;

};

#endif
