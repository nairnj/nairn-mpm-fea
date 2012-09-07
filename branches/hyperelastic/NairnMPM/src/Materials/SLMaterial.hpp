/********************************************************************************
    SLMaterial.hpp
    NairnMPM
    
    Created by John Nairn, 11/12/2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
	
	Add rate dependence to yield stress of the MGSCGLMaterial
	See D. J. Steinberg and C. M. Lund, J. Appl. Phys., v64, 1528-1533 (1989)
	"A constitutive model for strain rates from 10^-4 to 10^6 s^-1

	Dependencies
		MGSCGLMaterial.hpp Isoplasticity.hpp MaterialBase.hpp
********************************************************************************/

#ifndef SLMATERIAL

#define SLMATERIAL 18
#define YT_HISTORY 1
#define EPDOT_HISTORY 2
#define MAX_ITERATIONS 10
#define PRECISION_FACTOR 100000.

#include "Materials/MGSCGLMaterial.hpp"

class SLMaterial : public MGSCGLMaterial
{
    public:
		// unique properties
		double Uk,YP,C1,C2;
        
        // constructors and destructors
		SLMaterial();
		SLMaterial(char *matName);
		
		// initialize
        virtual char *InputMat(char *,int &);
		virtual void InitialLoadMechProps(int,int);
		virtual void PrintYieldProperties(void);
		virtual char *MaterialData(void);
        virtual const char *VerifyProperties(int);
				
		// methods
		virtual double GetPressureChange(MPMBase *,double &,int);
		virtual double GetYield(MPMBase *,int,double);
 		virtual double GetKPrime(MPMBase *,int,double);
		double GetEpdot(double);
		virtual double SolveForLambdaBracketed(MPMBase *,int,double,Tensor *,double);
		virtual void ElasticUpdateFinished(MPMBase *,int,double);
				
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
		virtual double GetHistory(int,char *);
 		
    protected:
		double YPred,TwoUkkT,C2red,currentYTred,constantYT;
		double epdotmin,epdotmax,YTmin,YTprecision;
		bool isConstantYT;

};

#endif

