/********************************************************************************
    NewMaterial.hpp
    NairnMPM
    
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
		//virtual bool ValidateForUse(int) const;
        //virtual void SetInitialParticleState(MPMBase *,int) const;
	
		// History-dependent properties
		//virtual char *InitHistoryData(void);
		//virtual double GetHistory(int,char *) const;
		
		// step methods
		//virtual voidv GetTransportProps(MPMBase *,int,TransportProperties *) const;
		//virtual double GetHeatCapacity(MPMBase *) const;
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
				
		// accessors
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag() const;
		virtual double WaveSpeed(bool,MPMBase *) const;
		//virtual double ShearWaveSpeed(bool,MPMBase *) const;
        //virtual double CurrentWaveSpeed(bool,MPMBase *) const;
		//virtual double MaximumDiffusion(void) const;
        //virtual double MaximumDiffusivity(void) const;
        //Tensor NewMaterial::GetStress(Tensor *sp,double pressure) const;
        //bool PartitionsElasticAndPlasticStrain(void);
        //bool SupportsArtificialViscosity(void) const;
		
    protected:

};

#endif

