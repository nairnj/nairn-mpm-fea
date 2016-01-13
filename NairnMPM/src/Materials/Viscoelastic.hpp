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

// This option switches to using Kirchoff stress even in this small strain material
// The shear parts are not implemented yet
#define USE_KIRCHOFF_STRESS

enum {XX_HISTORY=0,YY_HISTORY,XY_HISTORY,ZZ_HISTORY,XZ_HISTORY,YZ_HISTORY};
enum {MGJ_HISTORY,MGJRES_HISTORY};
enum {LINEAR_PRESSURE=0,MGEOS_PRESSURE};

class Viscoelastic : public MaterialBase
{
    public:
	
        // constructors and destructors
        Viscoelastic();
        Viscoelastic(char *);
        
        // innitialization methods
        virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
	
		// history data
		virtual char *InitHistoryData(char *,MPMBase *);
		virtual double GetHistory(int,char *) const;
	
		// const methods
		virtual void PrintMechanicalProperties(void) const;
		virtual void ValidateForUse(int) const;
    
		// methods
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int) const;
        virtual double GetCpMinusCv(MPMBase *) const;
#ifdef USE_KIRCHOFF_STRESS
		virtual void UpdatePressure(MPMBase *,double,ResidualStrains *,double,double,double,double,double &,double &) const;
#else
		virtual void UpdatePressure(MPMBase *,double,ResidualStrains *,double,const Matrix3 *,double,double &,double &) const;
#endif
		
		// accessors
        virtual Vector ConvertJToK(Vector,Vector,Vector,int);
        virtual double WaveSpeed(bool,MPMBase *) const;
		virtual const char *MaterialType() const;
		virtual int MaterialTag() const;
		virtual Tensor GetStress(Tensor *,double,MPMBase *) const;
		virtual bool SupportsArtificialViscosity(void) const;
		virtual double CurrentWaveSpeed(bool,MPMBase *,int) const;
#ifdef USE_KIRCHOFF_STRESS
		virtual double GetCurrentRelativeVolume(MPMBase *,int) const;
#endif
		
    private:
		int pressureLaw;
		double G0,K,aI;
		// double betaI;        // defined in superclass
		int ntaus;
		double *Gk,*tauk;
	
        int currentGk,currentTauk;
		double CTE,CME,Ka2sp;
		double Kered,Gered,*TwoGkred;
        
		// MG-EOS properties
		double gamma0,C0,S1,S2,S3;
		double C0squared;
		int mptrHistory;
};

#endif

