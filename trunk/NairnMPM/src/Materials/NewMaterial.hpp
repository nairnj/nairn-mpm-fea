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
        virtual const char *VerifyProperties(int);
		//virtual void InitialLoadMechProps(int,int);
		//virtual void InitialLoadTransProps(void);
		virtual void PrintMechanicalProperties(void);
		//virtual void PrintTransportProperties(void);
		//virtual bool ValidateForUse(int);
	
		// History-dependent properties
		//virtual char *MaterialData(void);
		//virtual double GetHistory(int,char *);
		
		// step methods
        //virtual void LoadMechanicalProps(MPMBase *,int);
		//virtual void LoadTransportProps(MPMBase *,int);
		//virtual double GetHeatCapacity(MPMBase *);
		//virtual double GetHeatCapacityVol(MPMBase *);
		virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,int);
		virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int);
				
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
		virtual double WaveSpeed(bool,MPMBase *);
		//virtual double ShearWaveSpeed(bool,MPMBase *);
        //virtual double CurrentWaveSpeed(bool,MPMBase *);
		//virtual double MaximumDiffusion(void);
        //virtual double MaximumDiffusivity(void);
		
    protected:

};

#endif

