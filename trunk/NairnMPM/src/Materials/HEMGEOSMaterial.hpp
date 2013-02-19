/********************************************************************************
	MGSCGLMaterial.hpp
	NairnMPM

	Created by John Nairn, 11/12/2008.
	Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		HEIsotropic.hpp (HyperElastic.hpp, MaterialBase.hpp)
********************************************************************************/

#ifndef HEMGEOSMATERIAL

#define HEMGEOSMATERIAL 25

#include "Materials/HEIsotropic.hpp"

class HEMGEOSMaterial : public HEIsotropic
{
	public:
		// unique properties
		double gamma0,C0,S1,S2,S3;
		
		// constructors and destructors
		HEMGEOSMaterial();
		HEMGEOSMaterial(char *matName);
		
		// initialize
		virtual char *InputMat(char *,int &);
		virtual const char *VerifyProperties(int);
		virtual void InitialLoadMechProps(int,int);
		virtual void PrintMechanicalProperties(void);
		virtual void PrintTransportProperties(void);
		virtual void ValidateForUse(int);
		
		// methods
		virtual void UpdatePressure(MPMBase *,double,double,int,double,double);
		virtual double GetCurrentRelativeVolume(MPMBase *);
		
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
		virtual double CurrentWaveSpeed(bool,MPMBase *);
		virtual bool SupportsArtificialViscosity(void);
		
	protected:
		double Keffred,Gratio;
		double C0squared;
	
};

#endif
