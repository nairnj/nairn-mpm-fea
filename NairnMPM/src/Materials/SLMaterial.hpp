/********************************************************************************
    SLMaterial.hpp
    nairn-mpm-fea
    
    Created by John Nairn, 11/12/2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
	
	Add rate dependence to yield stress
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

// plastic law properties
typedef struct {
	double Gratio;
	double TwoUkkT;
	double currentYTred;
	double constantYT;
	bool isConstantYT;
	double epdotmin;
	double epdotmax;
} SLProperties;

class SLMaterial : public SCGLHardening
{
    public:
        
        // constructors and destructors
		SLMaterial();
		SLMaterial(MaterialBase *);
		
		// initialize
        virtual char *InputHardeningProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
        virtual void PrintYieldProperties(void) const;
	
		// history data
		virtual int HistoryDoublesNeeded(void) const;
		
		// copy of properties
        virtual int SizeOfHardeningProps(void) const;
		virtual void *GetCopyOfHardeningProps(MPMBase *,int,void *,int);
		virtual void DeleteCopyOfHardeningProps(void *,int) const;
		virtual double GetShearRatio(MPMBase *,double,double,void *,int) const;
	
		// methods
        virtual double GetYield(MPMBase *,int,double,HardeningAlpha *,void *) const;
        virtual double GetKPrime(MPMBase *,int,double,HardeningAlpha *,void *) const;
        virtual double GetK2Prime(MPMBase *,double,double,HardeningAlpha *a,void *) const;
	
		// return mapping
        virtual double SolveForLambdaBracketed(MPMBase *,int,double,Tensor *,double,
											   double,double,double,HardeningAlpha *a,void *,int) const;
		double GetEpdot(double,double) const;
   
		// update
		virtual void ElasticUpdateFinished(MPMBase *,int,double,int) const;
	
		// accessors
        virtual const char *GetHardeningLawName(void) const;
 		
    protected:
		// unique properties
		double UkOverk,YP,C1,C2;
	
		// independent of particle state
		double YPred,C2red,YTmin,YTprecision;
	
};

#endif

