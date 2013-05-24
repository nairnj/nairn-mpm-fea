/********************************************************************************
    HEAnisotropic.hpp
    nairn-mpm-fea
    
    Created by John Nairn, Sept 27, 2011.
    Copyright (c) 2011 John A. Nairn, All rights reserved.

	Dependencies
		HyperElastic.hpp (MaterialBase.hpp)
********************************************************************************/

#ifndef HEANISOTROPIC

#define HEANISOTROPIC 21

#define J_HISTORY 0
#define THETA12_HISTORY 1

#include "Materials/HyperElastic.hpp"

class HEAnisotropic : public HyperElastic
{
    public:
        // constructors and destructors
		HEAnisotropic();
		HEAnisotropic(char *matName);
		
		// initialize
        virtual char *InputMat(char *,int &);
        virtual const char *VerifyAndLoadProperties(int);
        void ValidateForUse(int np) const;
		virtual void PrintMechanicalProperties(void) const;
        char *InitHistoryData(void);
		
		// step methods
        void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
				
		// accessors
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag() const;
		virtual double WaveSpeed(bool,MPMBase *) const;
		
    protected:
		// unique properties
	
};

#endif

