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

enum { AP_LINEAR=0,AP_NONLINEAR1,AP_NONLINEAR2,AP_EXPONENTIAL,AP_UNKNOWN };
#define ALPHA_EPS 1.e-12

// plastic law properties
typedef struct {
	double aint;
	double sAsmag;
	Tensor dfds;
	double dfCdf;
	Tensor Cdf;
	ElasticProperties *ep;
} AnisoPlasticProperties;

// Hill properties for use with reduced stresses
typedef struct {
    double F;
    double G;
    double H;
    double L;       // actually 2L
    double M;       // actually 2L
    double N;       // actually 2L
    double syxx2;   // 1/sig(xx)^2
    double syyy2;   // 1/sig(yy)^2
    double syzz2;   // 1/sig(zz)^2
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
        static void GetHillDfDsigma(Tensor &stk,int np,AnisoPlasticProperties *p,const HillProperties *h);
    
   protected:
		double syxx,syyy,syzz,tyyz,tyxz,tyxy;
		double sqrt23OversigmaYref;
        HillProperties h;

};

#endif
