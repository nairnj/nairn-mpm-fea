/********************************************************************************
    NewMaterial.hpp
    nairn-mpm-fea
    
    Created by John Nairn, July 13, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef NEWMATERIAL

#define NEWMATERIAL 100

#include "Materials/MaterialBase.hpp"

class NewMaterial : public MaterialBase
{
    public:
		// unique properties
		double newproperty;
        
        // constructors and destructors
		NewMaterial();
		NewMaterial(char *matName);
		
		// initialize
        virtual char *InputMat(char *,int &);
        virtual const char *VerifyAndLoadProperties(int);
		//virtual void FillTransportProperties(TransportProperties *);
		virtual void PrintMechanicalProperties(void) const;
		//virtual void PrintTransportProperties(void) const;
		//virtual void ValidateForUse(int) const;
	
		// History-dependent properties
        //virtual void SetInitialParticleState(MPMBase *,int) const;
		//virtual char *InitHistoryData(void);
		//virtual double GetHistory(int,char *) const;
		
		// contitutive law methods
        //virtual int SizeOfMechanicalProperties(int &altBufferSize) const
        //virtual void *MaterialBase::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer) const
		//virtual void GetTransportProps(MPMBase *,int,TransportProperties *) const;
		//virtual double GetHeatCapacity(MPMBase *) const;
		//virtual double NewMaterial::GetCpMinusCv(MPMBase *) const;
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
				
		// accessors
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag() const;
		virtual double WaveSpeed(bool,MPMBase *) const;
		//virtual double ShearWaveSpeed(bool,MPMBase *) const;
        //virtual double CurrentWaveSpeed(bool,MPMBase *) const;
		//virtual double MaximumDiffusion(void) const;
        //virtual double MaximumDiffusivity(void) const;
        //virtual Tensor GetStress(Tensor *sp,double pressure) const;
        //virtual bool PartitionsElasticAndPlasticStrain(void);
        //virtual bool SupportsArtificialViscosity(void) const;
		
    protected:

};

#endif

