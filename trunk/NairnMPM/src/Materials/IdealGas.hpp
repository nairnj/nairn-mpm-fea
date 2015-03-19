/********************************************************************************
    IdealGas.hpp
    nairn-mpm-fea
    
    Created by Edward Le and Peter Mackenzie on May 26, 2012.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
 
	Ideal gas as a hyperelastic material. 
 
	Dependencies
		HyperElastic.hpp
********************************************************************************/

#ifndef IDEALGASMATERIAL

#define IDEALGASMATERIAL 22

#include "Materials/HyperElastic.hpp"

enum {P0_PROP=0,RHO0_PROP,T0_PROP,CTE_PROPS,GAS_PROPS};

class IdealGas : public HyperElastic
{
    public:
 	    double P0;    // initial pressure
	    double T0;    // initial temperature in Kelvin
	
        // constructors and destructors
        IdealGas();
        IdealGas(char *);
        
        // initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
	
		// const methods
		virtual void ValidateForUse(int) const;
		virtual void PrintMechanicalProperties(void) const;
        void SetInitialParticleState(MPMBase *,int) const;
 		
		// methods
		virtual double GetCpMinusCv(MPMBase *) const;
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
    
		// accessors
		virtual double WaveSpeed(bool,MPMBase *) const;
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag() const;
		
    private:
	    double P0sp;  // mass specific initial pressure
		double CpMinusCv;
		double gammaAdiabatic;
};

#endif
