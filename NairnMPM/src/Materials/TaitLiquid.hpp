/********************************************************************************
    TaitLiquid.hpp
    nairn-mpm-fea

    Created by John Nairn, Dec 4, 2013.
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Dependencies
        HyperElastic.hpp,MaterialBase.hpp
 ********************************************************************************/

#ifndef TAITLIQUID

#define TAITLIQUID 27
#define TAIT_C 0.0894

class Expression;

#include "Materials/HyperElastic.hpp"

class TaitLiquid : public HyperElastic
{
    public:
    
        // constructors and destructors
        TaitLiquid(char *,int);
        
        // initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
        virtual void PrintMechanicalProperties(void) const;
        virtual void ValidateForUse(int) const;
        virtual void SetInitialParticleState(MPMBase *,int,int) const;
    
        // history data
		virtual int SizeOfHistoryData(void) const;
 		virtual char *InitHistoryData(char *,MPMBase *);
   		virtual int NumberOfHistoryDoubles(void) const;
    
        // contitutive law methods
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int) const;
		virtual double GetTwoEtaOverRho(double) const;
		virtual double BracketContactLawShearRate(double,double,double &,double &,double &,double &) const;
    
        // accessors
        virtual const char *MaterialType(void) const;
        virtual double WaveSpeed(bool,MPMBase *) const;
        virtual double CurrentWaveSpeed(bool,MPMBase *,int) const;
        virtual Tensor GetStress(Tensor *sp,double pressure,MPMBase *) const;
		virtual void SetStress(Tensor *,MPMBase *) const;
		virtual void IncrementThicknessStress(double dszz,MPMBase *mptr) const;
		virtual double GetCurrentRelativeVolume(MPMBase *,int) const;
		virtual void SetPressureFunction(char *);
		virtual double GetViscosity(double shearRate) const;
		virtual bool SupportsArtificialViscosity(void) const;
	
    protected:
		// unique properties
		vector<double> viscosity;
		vector<double> logShearRate;
		int numViscosity;
	
        double *TwoEtasp;
		double gamma0;
		Expression *function;
    
};

#endif

