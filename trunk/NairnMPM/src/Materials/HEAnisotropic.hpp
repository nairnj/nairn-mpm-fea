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
    double G1,G2,E,E1,E2,G3,G1sp,G2sp,G3sp,THETA0,E1sp,E2sp;
    
        // constructors and destructors
		HEAnisotropic();
		HEAnisotropic(char *matName);
		
		// initialize
        virtual char *InputMat(char *,int &);
        virtual const char *VerifyProperties(int);
		virtual void PrintMechanicalProperties(void);
		virtual void InitialLoadMechProps(int,int);
        virtual char *InitHistoryData(void);
		
		// step methods
    virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int);
				
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
		virtual double WaveSpeed(bool,MPMBase *);
		
    protected:
		// unique properties
	
};

#endif

