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
		virtual char *InputMaterialProperty(char *,int &,double &);
		char *InitHistoryData(void);
		
		// const methods
		virtual void PrintYieldProperties(void) const;
	
		// hardening law terms
		virtual double GetDfAlphaDotH(MPMBase *,int,AnisoHardProperties *p) const;
		virtual void UpdatePlasticInternal(MPMBase *,int,AnisoHardProperties *p) const;
		virtual void UpdateTrialAlpha(MPMBase *,int,AnisoHardProperties *) const;
		virtual void UpdateTrialAlpha(MPMBase *,int,double,AnisoHardProperties *) const;
		virtual double GetYield(AnisoHardProperties *p) const;
	
		// accessors
		virtual const char *MaterialType(void) const;
		double GetHistory(int,char *) const;
		virtual int MaterialTag() const;
};

#endif
