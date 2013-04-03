/********************************************************************************
    IdealGas.hpp
    NairnMPM
    
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
        virtual char *InputMat(char *,int &);
		virtual const char *VerifyProperties(int);
		virtual void ValidateForUse(int);
        virtual void InitialLoadMechProps(int,int);
		virtual void PrintMechanicalProperties(void);
        void SetInitialParticleState(MPMBase *,int);
 		
		// methods
        virtual void MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np);
    
		// accessors
		virtual double WaveSpeed(bool,MPMBase *);
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
        virtual double CurrentWaveSpeed(bool,MPMBase *);
		
    private:
	    double P0sp;  // mass specific initial pressure
	    double Psp;   // current specific pressure
	    double Temp;  // current temperature (in Kelvin)
};

#endif
