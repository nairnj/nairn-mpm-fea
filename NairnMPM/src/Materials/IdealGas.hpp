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
	
        // constructors and destructors
        IdealGas(char *,int);
        
        // initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
	
		// const methods
		virtual void ValidateForUse(int) const;
		virtual void PrintMechanicalProperties(void) const;
        void SetInitialParticleState(MPMBase *,int,int) const;
 		
		// history data
		virtual int SizeOfHistoryData(void) const;
		virtual char *InitHistoryData(char *,MPMBase *);
		virtual void ResetHistoryData(char *,MPMBase *);
		virtual int NumberOfHistoryDoubles(void) const;

		// methods
		virtual double GetCpMinusCv(MPMBase *) const;
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int) const;

		// accessors
        virtual Tensor GetStress(Tensor *sp,double pressure,MPMBase *) const;
        virtual void SetStress(Tensor *,MPMBase *) const;
        virtual void IncrementThicknessStress(double dszz,MPMBase *mptr) const;
		virtual double WaveSpeed(bool,MPMBase *) const;
        virtual double CurrentWaveSpeed(bool,MPMBase *,int) const;
		virtual const char *MaterialType(void) const;
		virtual bool SupportsDiffusion(void) const;
		virtual bool SupportsArtificialViscosity(void) const;
	
    private:
		double P0;    // initial pressure
		double T0;    // initial temperature in Kelvin
        double vdwa;    // van der Waals a terms
        double vdwb;    // van der Waals a terms
        double mu0;     // reference viscosity
        double cSuth;       // Sutherlands constant

	    double P0sp;  // mass specific initial pressure
		double CpMinusCv;
		double gammaAdiabatic;
        bool vanderWaalsGas;
        double aprime,bprime,V0n,P0prime;
        double twoMu0sp;
        bool realPressure;
};

#endif
