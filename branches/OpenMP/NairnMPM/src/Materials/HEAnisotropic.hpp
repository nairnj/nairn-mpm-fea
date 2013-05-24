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
		virtual void PrintMechanicalProperties(void) const;
		
		// step methods
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,int,void *,ResidualStrains *) const;
		virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int,void *,ResidualStrains *) const;
				
		// accessors
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag() const;
		virtual double WaveSpeed(bool,MPMBase *) const;
		
    protected:
		// unique properties
	
};

#endif

