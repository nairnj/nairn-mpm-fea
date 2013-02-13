/********************************************************************************
    AnisoPlasticity.hpp
    NairnMPM
    
    Created by John Nairn, June 12, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		Orthotropic.hpp (TranIsotropic.hpp, Elastic.hpp, MaterialBase.hpp)
********************************************************************************/

#ifndef _ANISOPLASTICITY_

#define _ANISOPLASTICITY_

#include "Materials/Orthotropic.hpp"

class AnisoPlasticity : public Orthotropic
{
    public:
        
        // constructors and destructors
		AnisoPlasticity();
		AnisoPlasticity(char *matName);
		
        // initialize
        virtual char *InputMat(char *,int &);
		virtual const char *VerifyProperties(int);
		virtual void ValidateForUse(int);
		virtual void InitialLoadMechProps(int,int);
		virtual void PrintMechanicalProperties(void);
		virtual void PrintYieldProperties(void);
			
		// methods
		virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,int);
		virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int);
		virtual double SolveForLambdaAP(MPMBase *mptr,int,double,Tensor *);
		virtual double GetFkFromLambdak(MPMBase *,Tensor *,Tensor *,double,int);
		virtual void GetDfCdf(MPMBase *,Tensor *,int);
		virtual void UpdateStress(Tensor *,Tensor *,double,int);
 		
		// subclass must provide platic potential functions
		virtual void UpdateTrialAlpha(MPMBase *,int) = 0;
		virtual void UpdateTrialAlpha(MPMBase *,int,double) = 0;
		virtual double GetF(MPMBase *,Tensor *,int) = 0;
		virtual void GetDfDsigma(MPMBase *,Tensor *,int) = 0;
		virtual double GetDfAlphaDotH(MPMBase *,int,Tensor *) = 0;
		virtual void UpdatePlasticInternal(MPMBase *,int) = 0;
    
        // accessors
        virtual bool PartitionsElasticAndPlasticStrain(void);
		
   protected:
		double syxx,syyy,syzz,tyyz,tyxz,tyxy;
		double syxxred2,syyyred2,syzzred2,tyyzred2,tyxzred2,tyxyred2;		// equal to 1/yield^2 and reduced
		double dfdsxx,dfdsyy,dfdszz,dfdtyz,dfdtxz,dfdtxy;
		double dfdsxxrot,dfdsyyrot,dfdszzrot,dfdtyzrot,dfdtxzrot,dfdtxyrot;
		double Cdfxx,Cdfyy,Cdfzz,Cdfyz,Cdfxz,Cdfxy,dfCdf;
		double cos2t,sin2t,costsint;										// 2D rotation terms calculated once per step
		double rzyx[6][6];													// 3D rotation matrix calcualted once per step

};

#endif
