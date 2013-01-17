/********************************************************************************
    MGSCGLMaterial.hpp
    NairnMPM
    
    Created by John Nairn, 11/12/2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		Isoplasticity.hpp MaterialBase.hpp
********************************************************************************/

#ifndef MGSCGLMATERIAL

#define MGSCGLMATERIAL 17

#include "Materials/IsoPlasticity.hpp"

class MGSCGLMaterial : public IsoPlasticity
{
    public:
		// unique properties
		double gamma0,GPp,GTp,beta,nhard,yieldMax;
		double C0,S1,S2,S3;
        
        // constructors and destructors
		MGSCGLMaterial();
		MGSCGLMaterial(char *matName);
		
		// initialize
        virtual char *InputMat(char *,int &);
		virtual const char *VerifyProperties(int);
		virtual void InitialLoadMechProps(int,int);
		virtual void PrintMechanicalProperties(void);
		virtual void PrintYieldProperties(void);
		virtual void PrintTransportProperties(void);
		virtual void ValidateForUse(int);
	
		// methods
        void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int);
		virtual double GetPressureChange(MPMBase *,double &,double,int);
		virtual double GetYield(MPMBase *,int,double);
 		virtual double GetKPrime(MPMBase *,int,double);
		virtual double GetK2Prime(MPMBase *,double,double);
        virtual double GetCurrentRelativeVolume(MPMBase *);
				
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
        virtual double WaveSpeed(bool,MPMBase *);
        virtual double CurrentWaveSpeed(bool,MPMBase *);
		
    protected:
        double GPpred,G0red,Keffred,Gratio;
        double yldMaxred;
        double C0squared;

};

#endif

