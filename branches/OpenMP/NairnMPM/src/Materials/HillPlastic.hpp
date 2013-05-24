/********************************************************************************
    HillPlastic.hpp
    nairn-mpm-fea
    
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
		char *InitHistoryData(void);
		
		// const methods
		virtual void PrintYieldProperties(void) const;
	
		// hardening law terms
		virtual double GetDfAlphaDotH(MPMBase *,int,Tensor *,AnisoPlasticProperties *p) const;
		virtual void UpdatePlasticInternal(MPMBase *,int,AnisoPlasticProperties *p) const;
		virtual void UpdateTrialAlpha(MPMBase *,int,AnisoPlasticProperties *) const;
		virtual void UpdateTrialAlpha(MPMBase *,int,double,AnisoPlasticProperties *) const;
		virtual double GetYield(AnisoPlasticProperties *p) const;
	
		// accessors
		virtual const char *MaterialType(void) const;
		double GetHistory(int,char *) const;
		virtual int MaterialTag() const;
};

#endif
