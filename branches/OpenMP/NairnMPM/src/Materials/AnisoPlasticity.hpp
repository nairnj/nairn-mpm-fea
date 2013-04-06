/********************************************************************************
    AnisoPlasticity.hpp
    NairnMPM
    
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
	ElasticProperties *ep;
	double aint;
	double minush;
	Tensor dfds;
	double dfCdf;
	Tensor Cdf;
	double rzyx[6][6];			// 3D rotation matrix calcualted once per step
} AnisoPlasticProperties;

class AnisoPlasticity : public Orthotropic
{
    public:
        
        // constructors and destructors
		AnisoPlasticity();
		AnisoPlasticity(char *matName);
		
        // initialize
        virtual char *InputMat(char *,int &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void ValidateForUse(int) const;
		virtual void PrintMechanicalProperties(void) const;
		virtual void PrintYieldProperties(void) const;
			
		// methods
		void *GetCopyOfMechanicalProps(MPMBase *,int) const;
		void DeleteCopyOfMechanicalProps(void *,int) const;
		virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,int,void *,ResidualStrains *) const;
		virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int,void *,ResidualStrains *) const;
		virtual double SolveForLambdaAP(MPMBase *mptr,int,double,Tensor *,AnisoPlasticProperties *p) const;
		virtual double GetFkFromLambdak(MPMBase *,Tensor *,Tensor *,double,int,AnisoPlasticProperties *) const;
		virtual void UpdateStress(Tensor *,Tensor *,double,int,AnisoPlasticProperties *p) const;
	
		// Hill methods
		virtual double GetMagnitudeRotatedHill(Tensor *,Tensor *,int,AnisoPlasticProperties *) const;
		virtual void GetDfCdf(Tensor *,int,AnisoPlasticProperties *p) const;
		virtual void GetDfDsigma(Tensor *,int,AnisoPlasticProperties *p) const;
 		
		// hardening term methods (move to hardening law class when want more hardening options)
		virtual void UpdateTrialAlpha(MPMBase *,int,AnisoPlasticProperties *) const = 0;
		virtual void UpdateTrialAlpha(MPMBase *,int,double,AnisoPlasticProperties *p) const = 0;
		virtual double GetYield(AnisoPlasticProperties *p) const = 0;
		virtual double GetDfAlphaDotH(MPMBase *,int,Tensor *,AnisoPlasticProperties *p) const = 0;
		virtual void UpdatePlasticInternal(MPMBase *,int,AnisoPlasticProperties *p) const = 0;
    
        // accessors
        virtual bool PartitionsElasticAndPlasticStrain(void);
		
   protected:
		double syxx,syyy,syzz,tyyz,tyxz,tyxy;
		double syxxred2,syyyred2,syzzred2,tyyzred2,tyxzred2,tyxyred2;		// equal to 1/yield^2 and reduced
		double fTerm,gTerm,hTerm;

};

#endif
