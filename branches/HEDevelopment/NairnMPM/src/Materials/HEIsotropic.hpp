/********************************************************************************
    HEIsotropic.hpp
    NairnMPM
    
    Created by John Nairn, Sept 27, 2011.
    Copyright (c) 2011 John A. Nairn, All rights reserved.

	Dependencies
		HyperElastic.hpp (MaterialBase.hpp)
********************************************************************************/

#ifndef HEISOTROPIC

#define HEISOTROPIC 24

#include "Materials/HyperElastic.hpp"

enum {G11_PROP=0,KBULKK_PROP,CTEE_PROP,HEISOTROPIC_PROPS};

class HEIsotropic : public HyperElastic
{
    public:
    double G1,Kbulk;
        // constructors and destructors
		HEIsotropic();
		HEIsotropic(char *matName);
		
		// initialize
        virtual char *InputMat(char *,int &);
        virtual const char *VerifyProperties(int);
		virtual void PrintMechanicalProperties(void);
		virtual void InitialLoadMechProps(int,int);
        virtual char *MaterialData(void);
        virtual void SetInitialParticleState(MPMBase *,int);
		
		// step methods
		virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,int);
		virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int);
        virtual double GetVolumetricTerms(double,double *);
        
				
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
		virtual double WaveSpeed(bool,MPMBase *);
		
    protected:
		// unique properties
private:
    double G1sp, Ksp;
	
};

#endif

