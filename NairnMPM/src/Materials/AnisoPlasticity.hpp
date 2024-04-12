/********************************************************************************
    AnisoPlasticity.hpp
    nairn-mpm-fea
    
    Created by John Nairn, June 12, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		Orthotropic.hpp (TranIsotropic.hpp, Elastic.hpp, MaterialBase.hpp)
********************************************************************************/

#ifndef _ANISOPLASTICITY_

#define _ANISOPLASTICITY_

#include "Materials/Orthotropic.hpp"

// Hill styles
enum { SQRT_TERMS=1,SQUARED_TERMS };

enum { AP_LINEAR=0,AP_NONLINEAR1,AP_NONLINEAR2,AP_EXPONENTIAL,AP_UNKNOWN };
#define ALPHA_EPS 1.e-12

// plastic law properties
typedef struct {
	double aint;
	Tensor dfdsPsigma;			// dfds or P sigma
	Tensor CdfQsigma;			// C dfds or Q sigma
	double sAQsmag;				// sqrt(sAs) or sqrt(sQs)
	double dfCdf;				// only for square root terms
	ElasticProperties *ep;
} AnisoPlasticProperties;

// Hill properties for use with reduced stresses
typedef struct {
    double F;
    double G;
    double H;
    double L;       // actually 2L
    double M;       // actually 2M
    double N;       // actually 2N
    double syxx2;   // 1/sig(xx)^2 = G+H
    double syyy2;   // 1/sig(yy)^2 = F+H
    double syzz2;   // 1/sig(zz)^2 = F+G
	// NEW_HILL terms
	int style;		// SQRT_TERMS (1) ot SQUARED_TERMS (2)
	double CP11;
	double CP12;
	double CP13;
	double CP21;
	double CP22;
	double CP23;
	double CP31;
	double CP32;
	double CP33;
	double CP44;
	double CP55;
	double CP66;
	double Q11;
	double Q12;
	double Q13;
	double Q22;
	double Q23;
	double Q33;
	double Q44;
	double Q55;
	double Q66;
} HillProperties;

class AnisoPlasticity : public Orthotropic
{
    public:
		static int warnNonconvergence;
    
        // constructors and destructors
		AnisoPlasticity(char *,int);
		
        // initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		static void FillHillStyleProperties(int,HillProperties &,ElasticProperties &);
		virtual void ValidateForUse(int) const;
		virtual void PrintMechanicalProperties(void) const;
        virtual void PrintYieldProperties(void) const = 0;

		// methods
        virtual int SizeOfMechanicalProperties(int &) const;
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *,int) const;
		virtual int AltStrainContains(void) const;
		virtual ElasticProperties *GetElasticPropertiesPointer(void *) const;
	
		// Large Rotation Methods
		virtual void LRElasticConstitutiveLaw(MPMBase *,Matrix3 &,Matrix3 &,Matrix3 &,Matrix3 &,Matrix3 &,int,void *,ResidualStrains *) const;
		virtual Tensor SolveForPlasticIncrement(MPMBase *,int,double,Tensor &,AnisoPlasticProperties *) const;
		virtual void GetDfCdf(Tensor &,int,AnisoPlasticProperties *) const;
	
		// hardening term methods (move to hardening law class when want more hardening options)
		virtual double GetYield(AnisoPlasticProperties *p) const = 0;
		virtual double GetGPrime(AnisoPlasticProperties *) const = 0;
    
        // class function
        static void PrintAPYieldProperties(double,double,double,double,double,double);
        static double GetHillMagnitude(Tensor &,const HillProperties *,int);
        static void GetHillDfDsigmaPQsigma(Tensor &stk,int np,AnisoPlasticProperties *p,const HillProperties *h);
    
   protected:
		double syxx,syyy,syzz,tyyz,tyxz,tyxy;
		double sqrt23OversigmaYref;
        HillProperties h;

};

#endif
