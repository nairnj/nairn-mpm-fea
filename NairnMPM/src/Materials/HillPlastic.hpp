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
        // Harding constantsL Y/Yref = (1 + Khard alpha^nhard)
		double Khard,nhard,exphard,Kexp;
        int hardStyle;
        
        // constructors and destructors
		HillPlastic(char *,int);
		
        // methods
		virtual char *InputMaterialProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);

		// history data
		virtual char *InitHistoryData(char *,MPMBase *);
   		virtual int NumberOfHistoryDoubles(void) const;
		
		// const methods
		virtual void PrintYieldProperties(void) const;
	
		// hardening law terms
		virtual double GetYield(AnisoPlasticProperties *p) const;
		virtual double GetGPrime(AnisoPlasticProperties *) const;
	
		// accessors
		virtual const char *MaterialType(void) const;
};

#endif
