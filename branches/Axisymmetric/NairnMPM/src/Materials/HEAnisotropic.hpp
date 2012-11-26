/********************************************************************************
    HEAnisotropic.hpp
    NairnMPM
    
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
        virtual const char *VerifyProperties(int);
		virtual void PrintMechanicalProperties(void);
		virtual void InitialLoadMechProps(int,int);
		
		// step methods
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,int);
		virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int);
				
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
		virtual double WaveSpeed(bool,MPMBase *);
		
    protected:
		// unique properties
	
};

#endif

