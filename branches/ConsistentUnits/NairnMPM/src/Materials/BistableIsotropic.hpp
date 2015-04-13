/********************************************************************************
    BistableIsotropic.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Apr 11 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.

	Dependencies
		IsotropicMat.hpp
********************************************************************************/

#ifndef BISTABLEISO

#define BISTABLEISO 10

#include "Materials/IsotropicMat.hpp"

enum {K0_PROP=0,KD_PROP,G0_PROP,GD_PROP,A0_PROP,AD_PROP,
            TRANSITION_PROP,DVCRIT_PROP,DVOFF_PROP,DIFF0_PROP,
			DIFFD_PROP,KCOND0_PROP,KCONDD_PROP,B0_PROP,BD_PROP,BISTABLE_PROPS};

// lowest bit is dilated, second bit is sheared
#define INITIAL_STATE 0
#define DEFORMED_STATE 1

#define DILATION_RULE 1
#define DISTORTION_RULE 2
#define VONMISES_RULE 3

class BistableIsotropic : public IsotropicMat
{
    public:
        double K0,Kd;		// bulk moduli two possible states
        double G0,Gd;		// shear moduli two possible states
        double a0,ad;		// CTE two possible states
        double beta0,betad;		// CME two possible states
        double dVcrit;		// critical volumetric or deviatoric strain or stress for transition
        double dVii;		// volumetric strain offest when dilated
        int rule;			// transition rule ID
        bool reversible;	// is it reversible?
		double diff0,diffd;		// diffusion constants
		double kCond0,kCondd;	// conductivity constants
        
        // constructors and destructors
        BistableIsotropic();
        BistableIsotropic(char *);
		
		// initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual char *InitHistoryData(void);
		virtual const char *CurrentProperties(short,int);
	
		// const methods
		virtual void ValidateForUse(int) const;
		virtual void PrintMechanicalProperties(void) const;
		virtual void PrintTransportProperties(void) const;
        
        // override methods
#ifdef USE_PSEUDOHYPERELASTIC
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
#else
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,int,void *,ResidualStrains *) const;
#endif
        virtual void GetTransportProps(MPMBase *,int,TransportProperties *) const;
		
		// accessors
        virtual double WaveSpeed(bool,MPMBase *) const;
        virtual double MaximumDiffusion(void) const;
        virtual double MaximumDiffusivity(void) const;
		virtual const char *MaterialType(void) const;
		virtual double GetHistory(int,char *) const;
		virtual int MaterialTag() const;
    
    private:
        char readbs[BISTABLE_PROPS];
        short transState;	// current state for transport properties
		ElasticProperties pr2;
		TransportProperties tr2;
};

#endif
