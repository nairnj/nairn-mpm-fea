/********************************************************************************
    IsoPlasticity.hpp
    NairnMPM
    
    Created by John Nairn, June 16, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		IsotropicMat.hpp (Elastic.hpp, MaterialBase.hpp)
********************************************************************************/

#ifndef _ISOPLASTICITY_

#define _ISOPLASTICITY_
#define SQRT_TWOTHIRDS 0.8164965809277260
#define TWOTHIRDS 0.6666666666666667
#define ONETHIRD 0.3333333333333333
#define SQRT_EIGHT27THS 0.5443310539518174

#include "Materials/IsotropicMat.hpp"

class IsoPlasticity : public IsotropicMat
{
    public:
		// one yield stress for isotropic, plastic materials
		double yield;
        
        // constructors and destructors
		IsoPlasticity();
		IsoPlasticity(char *matName);
		
		// initialize
        virtual char *InputMat(char *,int &);
        virtual const char *VerifyProperties(int);
		virtual void InitialLoadMechProps(int,int);
		virtual char *MaterialData(void);
		virtual void PrintYieldProperties(void) = 0;							// subclass must provide
		
		// methods
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,int);
		virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int);
		
		// custom methods: Find yield function and solve for lambda
		virtual double GetPressureChange(MPMBase *,double &,int);
		virtual double GetMagnitudeS(Tensor *st,int);
		virtual void GetDfDsigma(double,Tensor *,int);
		virtual double SolveForLambda(MPMBase *,int,double,Tensor *,double);
        virtual double SolveForLambdaBracketed(MPMBase *,int,double,Tensor *,double);
        virtual void BracketSolution(MPMBase *,int,double,Tensor *,double,double *,double *);
		virtual bool LambdaConverged(int,double,double);
		virtual double GetYield(MPMBase *,int,double) = 0;				// subclass must provide
		virtual double GetKPrime(MPMBase *,int,double) = 0;				// subclass must provide
		virtual double GetK2Prime(MPMBase *,double,double) = 0;		    // subclass must provide
		virtual void ElasticUpdateFinished(MPMBase *,int,double);
		
		// Default internal variable as cumlative plastic strain
		virtual void UpdateTrialAlpha(MPMBase *,int);
		virtual void UpdateTrialAlpha(MPMBase *,int,double,double);
		virtual void UpdatePlasticInternal(MPMBase *,int);
		
		// accessors
		virtual double GetHistory(int,char *);
        virtual bool HasPlasticStrainForGradient(void);
		
    protected:
		double yldred,Gred,Kred;
		double psRed,psLr2G,psKred;
		int readYield;
		double dfdsxx,dfdsyy,dfdtxy,dfdszz,dfdtxz,dfdtyz;
		double alpint,dalpha;		// internal alpha variable and current increment

};

#endif

