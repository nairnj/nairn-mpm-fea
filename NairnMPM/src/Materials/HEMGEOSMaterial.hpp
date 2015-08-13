/********************************************************************************
	HEMGEOSMaterial.hpp
	nairn-mpm-fea

	Created by John Nairn, 11/12/2008.
	Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		HEIsotropic.hpp (HyperElastic.hpp, MaterialBase.hpp)
********************************************************************************/

#ifndef HEMGEOSMATERIAL

#define HEMGEOSMATERIAL 25

// old code reverts to new hyperelastic material
#define MGEOSMATERIAL 17

#include "Materials/HEIsotropic.hpp"

class HEMGEOSMaterial : public HEIsotropic
{
	public:
		
		// constructors and destructors
		HEMGEOSMaterial();
		HEMGEOSMaterial(char *matName);
		
		// initialize
		virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
	
		// const methods
		virtual void PrintMechanicalProperties(void) const;
		virtual void PrintTransportProperties(void) const;
		virtual void ValidateForUse(int) const;
		
		// methods
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *) const;
		virtual void UpdatePressure(MPMBase *,double,double,int,double,double,HEPlasticProperties *,ResidualStrains *,double) const;
		
		// accessors
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag() const;
		virtual double CurrentWaveSpeed(bool,MPMBase *) const;
		
	protected:
		// unique properties
		double gamma0,C0,S1,S2,S3;
		double C0squared;
	
};

#endif
