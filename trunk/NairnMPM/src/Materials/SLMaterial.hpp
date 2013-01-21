/********************************************************************************
    SLMaterial.hpp
    NairnMPM
    
    Created by John Nairn, 11/12/2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
	
	Add rate dependence to yield stress of the MGSCGLMaterial
	See D. J. Steinberg and C. M. Lund, J. Appl. Phys., v64, 1528-1533 (1989)
	"A constitutive model for strain rates from 10^-4 to 10^6 s^-1

	Dependencies
		SCGLHardening.hpp, HardeningLawBase.hpp, MaterialBase.hpp
********************************************************************************/

#ifndef _SLMATERIAL_

#define _SLMATERIAL_

#define YT_HISTORY 1
#define EPDOT_HISTORY 2
#define MAX_ITERATIONS 10
#define PRECISION_FACTOR 100000.

#include "Materials/SCGLHardening.hpp"

class SLMaterial : public SCGLHardening
{
    public:
		// unique properties
		double Uk,YP,C1,C2;
        
        // constructors and destructors
		SLMaterial();
		SLMaterial(MaterialBase *);
		
		// initialize
        virtual char *InputMat(char *,int &);
        virtual const char *VerifyProperties(int);
        virtual void PrintYieldProperties(void);
        virtual void InitialLoadMechProps(int,int);
		//virtual char *MaterialData(void);
 				
		// methods
        virtual double GetShearRatio(MPMBase *,double,double);
        virtual double GetYield(MPMBase *,int,double);
        virtual double GetKPrime(MPMBase *,int,double);
        virtual double GetK2Prime(MPMBase *,double,double);
        double GetEpdot(double YT);
        virtual double SolveForLambdaBracketed(MPMBase *,int,double,Tensor *,double,double,double);
   
		// accessors
        virtual const char *GetHardeningLawName(void);
 		
    protected:
		double YPred,TwoUkkT,C2red,currentYTred,constantYT;
		double epdotmin,epdotmax,YTmin,YTprecision;
		bool isConstantYT;

};

#endif

