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
		virtual void InitialLoadMechProps(int,int);
		virtual void PrintMechanicalProperties(void);
		virtual void PrintYieldProperties(void);
		virtual void ValidateUse(int);
		virtual bool ThreeDMaterial(void);
			
		// methods
		virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,int);
		virtual double SolveForLambda(MPMBase *mptr,int,double,Tensor *);
		virtual bool LambdaConverged(int,double,double);
 		
		// subclass must provide platic potential functions
		virtual void UpdateTrialAlpha(MPMBase *,int) = 0;
		virtual void UpdateTrialAlpha(MPMBase *,int,double) = 0;
		virtual double GetF(MPMBase *,double,double,double,double,int) = 0;
		virtual void GetDfDsigma(MPMBase *,Tensor *,int) = 0;
		virtual double GetDfAlphaDotH(MPMBase *,int,Tensor *) = 0;
		virtual void UpdatePlasticInternal(MPMBase *,int) = 0;
		
   protected:
		double syxx,syyy,syzz,tyxy;
		double syxxred2,syyyred2,syzzred2,tyxyred2;		// equal to 1/yield^2 and reduced
		double dfdsxx,dfdsyy,dfdtxy,dfdszz;
		double dfdsxxrot,dfdsyyrot,dfdtxyrot,dfdszzrot;
		double cos2t,sin2t,costsint;					// rotation terms calculated ones per step

};

#endif
