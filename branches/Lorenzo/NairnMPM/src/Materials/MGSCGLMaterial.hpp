/********************************************************************************
    MGSCGLMaterial.hpp
    NairnMPM
    
    Created by John Nairn, 11/12/2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		IsoPlasticity.hpp (IsotropicMat.hpp, Elastic.hpp, MaterialBase.hpp)
********************************************************************************/

#ifndef MGEOSMATERIAL

#define MGEOSMATERIAL 17

#include "Materials/IsoPlasticity.hpp"

class MGSCGLMaterial : public IsoPlasticity
{
    public:
		// unique properties
		double gamma0,C0,S1,S2,S3;
        
        // constructors and destructors
		MGSCGLMaterial();
		MGSCGLMaterial(char *matName);
		
		// initialize
        virtual char *InputMat(char *,int &);
		virtual const char *VerifyProperties(int);
		virtual void InitialLoadMechProps(int,int);
		virtual void PrintMechanicalProperties(void);
		virtual void PrintTransportProperties(void);
		virtual void ValidateForUse(int);
	
		// methods
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int);
		virtual void UpdatePressure(MPMBase *,double &,double,int);
        virtual double GetCurrentRelativeVolume(MPMBase *);
				
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
        virtual double WaveSpeed(bool,MPMBase *);
        virtual double CurrentWaveSpeed(bool,MPMBase *);
		virtual bool SupportsArtificialViscosity(void);
	
    protected:
        double GPpred,G0red,Keffred,Gratio;
        double C0squared,QAVred;

};

#endif

