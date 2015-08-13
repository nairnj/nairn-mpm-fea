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

typedef struct {
	double aint;
	double minush;
} AnisoHardProperties;

// plastic law properties
typedef struct {
	AnisoHardProperties hp;
	Tensor dfds;
	double dfCdf;
	Tensor Cdf;
	Tensor Cdf0;
	double snorm;
	ElasticProperties *ep;
} LRAnisoPlasticProperties;

// plastic law properties
typedef struct {
	AnisoHardProperties hp;
	Tensor dfds;
	double dfCdf;
	Tensor Cdf;
	Tensor Cdf0;
	double snorm;
	ElasticProperties ep;
	double rzyx[6][6];			// 3D rotation matrix calcualted once per step
} AnisoPlasticProperties;

class AnisoPlasticity : public Orthotropic
{
    public:
		static int warnNonconvergence;
    
        // constructors and destructors
		AnisoPlasticity();
		AnisoPlasticity(char *matName);
		
        // initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void ValidateForUse(int) const;
		virtual void PrintMechanicalProperties(void) const;
		virtual void PrintYieldProperties(void) const;
			
		// methods
        virtual int SizeOfMechanicalProperties(int &) const;
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *) const;
		virtual int AltStrainContains(void) const;
		virtual ElasticProperties *GetElasticPropertiesPointer(void *) const;
	
		// Lare Rotation Methods
		virtual void LRElasticConstitutiveLaw(MPMBase *,Matrix3,Matrix3,Matrix3,Matrix3,Matrix3 *,int,void *,ResidualStrains *) const;
		virtual double GetMagnitudeHill(Matrix3 &,int) const;
		virtual double LRSolveForLambdaAP(MPMBase *mptr,int,double,Matrix3 &,Matrix3 &,LRAnisoPlasticProperties *p) const;
		virtual void LRGetDfCdf(Matrix3 &,int,LRAnisoPlasticProperties *p) const;
		virtual void LRGetDfDsigma(Matrix3 &,int,LRAnisoPlasticProperties *p) const;
		virtual void LRUpdateStress(Matrix3 &,Matrix3 &,double,int,LRAnisoPlasticProperties *p) const;
		virtual double LRGetFkFromLambdak(MPMBase *,Matrix3 &,Matrix3 &,double,int,LRAnisoPlasticProperties *) const;
		virtual double LRPrintFk(MPMBase *,Matrix3 &,Matrix3 &,double,int,LRAnisoPlasticProperties *,double,double) const;

		// small rotation methods
		virtual void SRConstitutiveLaw2D(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
		virtual void SRConstitutiveLaw3D(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
		virtual double GetMagnitudeRotatedHill(Tensor *,Tensor *,int,AnisoPlasticProperties *) const;
		virtual double SolveForLambdaAP(MPMBase *mptr,int,double,Tensor *,AnisoPlasticProperties *p) const;
		virtual void GetDfCdf(Tensor *,int,AnisoPlasticProperties *p) const;
		virtual void GetDfDsigma(Tensor *,int,AnisoPlasticProperties *p) const;
		virtual void UpdateStress(Tensor *,Tensor *,double,int,AnisoPlasticProperties *p) const;
		virtual double GetFkFromLambdak(MPMBase *,Tensor *,Tensor *,double,int,AnisoPlasticProperties *) const;
		virtual double SRPrintFk(MPMBase *,Tensor *,Tensor *,double,int,AnisoPlasticProperties *,double,double) const;
 		
		// hardening term methods (move to hardening law class when want more hardening options)
		virtual void UpdateTrialAlpha(MPMBase *,int,AnisoHardProperties *) const = 0;
		virtual void UpdateTrialAlpha(MPMBase *,int,double,AnisoHardProperties *p) const = 0;
		virtual double GetYield(AnisoHardProperties *p) const = 0;
		virtual double GetDfAlphaDotH(MPMBase *,int,AnisoHardProperties *p) const = 0;
		virtual void UpdatePlasticInternal(MPMBase *,int,AnisoHardProperties *p) const = 0;
		
   protected:
		double syxx,syyy,syzz,tyyz,tyxz,tyxy;
		double syxxred2,syyyred2,syzzred2,tyyzred2,tyxzred2,tyxyred2;		// equal to 1/yield^2 and reduced
		double fTerm,gTerm,hTerm;

};

#endif
