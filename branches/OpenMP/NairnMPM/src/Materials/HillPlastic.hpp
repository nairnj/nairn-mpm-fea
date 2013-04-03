/********************************************************************************
    HillPlastic.hpp
    NairnMPM
    
    Created by John Nairn, June 12, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		AnisoPlasticity.hpp (Orthotropic.hpp TranIsotropic.hpp Elastic.hpp MaterialBase.hpp)
********************************************************************************/

#ifndef HILLPLASTIC

#define HILLPLASTIC 15

#include "Materials/AnisoPlasticity.hpp"

class HillPlastic : public AnisoPlasticity
{
    public:
        // Harding constantsL Y = Y0(1 + Khard alpha)^nhard
		double Khard,nhard;
        
        // constructors and destructors
		HillPlastic();
		HillPlastic(char *matName);
		
        // methods
		char *InputMat(char *,int &);
		virtual const char *VerifyAndLoadProperties(int);
		char *InitHistoryData(void);
		
		// const methods
		virtual void PrintYieldProperties(void) const;
	
		// plastic potential functions
		virtual void UpdateTrialAlpha(MPMBase *,int);
		virtual void UpdateTrialAlpha(MPMBase *,int,double);
		virtual double GetF(MPMBase *,Tensor *,int);
		virtual double GetF3D(MPMBase *,Tensor *);
		virtual void GetDfDsigma(MPMBase *,Tensor *,int);
		virtual void GetDfDsigma3D(MPMBase *,Tensor *);
		virtual double GetDfAlphaDotH(MPMBase *,int,Tensor *);
		virtual void UpdatePlasticInternal(MPMBase *,int);
		
		// accessors
		virtual const char *MaterialType(void) const;
		double GetHistory(int,char *) const;
		virtual int MaterialTag() const;
				
    protected:
		double fTerm,gTerm,hTerm;
		
		// most recent calculation
		double sxrot,syrot,txyrot;
		
		// internal variable while solving for lambda
		double aint,minush;
};

#endif
