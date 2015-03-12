/********************************************************************************
    MGSCGLMaterial.hpp
    nairn-mpm-fea
    
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
        virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;
		virtual void PrintTransportProperties(void) const;
		virtual void ValidateForUse(int) const;
	
		// methods
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *) const;
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
		virtual void UpdatePressure(MPMBase *,double &,double,int,PlasticProperties *,ResidualStrains *,double) const;
        virtual double GetCurrentRelativeVolume(MPMBase *) const;
				
		// accessors
        virtual Vector ConvertJToK(Vector,Vector,Vector,int);
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag() const;
        virtual double WaveSpeed(bool,MPMBase *) const;
        virtual double CurrentWaveSpeed(bool,MPMBase *) const;
		virtual bool SupportsArtificialViscosity(void) const;
	
    protected:
        double C0squared;

};

#endif

