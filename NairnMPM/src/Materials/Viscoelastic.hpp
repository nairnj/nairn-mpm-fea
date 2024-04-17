/********************************************************************************
    Viscoelastic.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef VISCOELASTIC

#define VISCOELASTIC 7

#include "Materials/MaterialBase.hpp"

enum {XX_HISTORY=0,YY_HISTORY,XY_HISTORY,ZZ_HISTORY,XZ_HISTORY,YZ_HISTORY};
enum {MGJ_HISTORY=0,MGJRES_HISTORY};
enum {LINEAR_PRESSURE=0,MGEOS_PRESSURE,TIME_DEPENDENT_PRESSURE};

class Viscoelastic : public MaterialBase
{
    public:
		static int warnExcessiveX;
	
        // constructors and destructors
        Viscoelastic(char *,int);
        
        // innitialization methods
        virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
	
		// history data
		virtual char *InitHistoryData(char *,MPMBase *);
		virtual void ResetHistoryData(char *,MPMBase *);
		virtual int NumberOfHistoryDoubles(void) const;
	
		// const methods
		virtual void PrintMechanicalProperties(void) const;
		virtual void ValidateForUse(int) const;
    
		// methods
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int) const;
        void GetAlphaArgs(double,double,double &,double &) const;
        double GetPSArg(double,double) const;
        virtual double GetCpMinusCv(MPMBase *) const;
		virtual void UpdatePressure(MPMBase *,double,ResidualStrains *,double,double,double,double,double &,double &) const;
	
		// accessors
        virtual Vector ConvertJToK(Vector,Vector,Vector,int);
        virtual double WaveSpeed(bool,MPMBase *) const;
		virtual const char *MaterialType() const;
		virtual Tensor GetStress(Tensor *,double,MPMBase *) const;
		virtual void SetStress(Tensor *,MPMBase *) const;
		virtual void IncrementThicknessStress(double dszz,MPMBase *mptr) const;
		virtual bool SupportsArtificialViscosity(void) const;
		virtual double CurrentWaveSpeed(bool,MPMBase *,int) const;
		virtual double GetCurrentRelativeVolume(MPMBase *,int) const;
		virtual bool SupportsDiffusion(void) const;
    
        // class methods
        static double GetVertialShift(MPMBase *,double,vector<double>,vector<double>,
                                      double,vector<double>,vector<double>);
        static double GetEffectiveIncrement(MPMBase *,ResidualStrains *,double,double,double,double,
                                            double,double,double,double,double);
		static double GetDtScale(double);

    private:
		int pressureLaw;
		double G0,K,aI;
		// double betaI;        // defined in superclass
		int ntaus;
		double *Gk,*tauk;
        int ntausK;
        double *Kk,*tauKk,*Kkred;
        int currentKk,currentTauKk;
	
        int currentGk,currentTauk;
		double CTE,CME,Ka2sp;
		double Kered,Gered,*TwoGkred;
        
		// MG-EOS properties
		double gamma0,C0,S1,S2,S3;
		double C0squared;
		double Kmax,Xmax;
    
        // WLF parameters
        double Tref,C1,C1base10,C2;
        double mref,Cm1,Cm2,Cm1base10,Cm2base10;

        // history variables
        int numJHistory,numHistory;
    
        // vertical shifting
        vector<double> bTemp;
        vector<double> bTValue;
        vector<double> bConc;
        vector<double> bCValue;
};

#endif

