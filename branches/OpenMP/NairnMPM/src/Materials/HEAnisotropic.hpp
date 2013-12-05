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
#define NUM_HISTORY 2

#include "Materials/HyperElastic.hpp"

class HEAnisotropic : public HyperElastic
{
    public:
        // constructors and destructors
        double THETA0, Eyarn, nu;
        HEAnisotropic();
		HEAnisotropic(char *matName);
		
		// initialize
        virtual char *InputMat(char *,int &);
        virtual const char *VerifyAndLoadProperties(int);
        virtual void ValidateForUse(int np) const;
		virtual void PrintMechanicalProperties(void) const;
        char *InitHistoryData(void);
		
		// step methods
        void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
				
		// accessors
        virtual double GetHistory(int,char *) const;
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag() const;
		virtual double WaveSpeed(bool,MPMBase *) const;
		
    protected:
		// unique properties
    bool balanced;
    double T10, T11, T12, T13, T14;
    double T20, T21, T22, T23, T24;
    double SH0, SH1, SH2, SH3, SH4, SH5, SH6, SH7, SH8, SH9, SH10, SH11, SH12;
	
};

#endif

