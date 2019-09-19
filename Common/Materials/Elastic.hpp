/********************************************************************************
    Elastic.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 24 2007.
    Copyright (c) 2007 John A. Nairn, All rights reserved.

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef _ELASTIC_

#define _ELASTIC_

#include "Materials/MaterialBase.hpp"

#ifdef MPM_CODE
// The full stiffness matrix in C
// alpha and beta are thermal and moisture expansion/poroelasticity coefficients
// although some elements are used for other things
typedef struct {
	double C[6][6];
	double alpha[8];
	double beta[6];
} ElasticProperties;

class SofeningLaw;

// softrning history variables
enum { SOFT_DAMAGE_STATE=0,DELTANORMAL,DELTASHEAR,DELTASHEAR2,DAMAGENORMAL,DAMAGESHEAR,DAMAGESHEAR2,
	NORMALDIR1,NORMALDIR2,NORMALDIR3,GCSCALING,RELATIVE_STRENGTH,RELATIVE_TOUGHNESS,
	SOFT_NUMBER_HISTORY };
enum { RECTANGULAR_SURFACE=0,ELLIPTICAL_SURFACE};

#endif

class Elastic : public MaterialBase
{
    public:
        
        // constructors and destructors
        Elastic(char *,int);
		
		// initialize
		virtual char *InputMaterialProperty(char *,int &,double &);
#ifdef MPM_CODE
		virtual void PrintCommonProperties(void) const;
#endif
	
		// methods
		void FillUnrotatedElasticProperties(ElasticProperties *,int);
#ifdef MPM_CODE
		virtual void HypoIncrementDeformation(MPMBase *,Matrix3) const;
        virtual double GetCpMinusCv(MPMBase *) const;
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int) const;
		virtual void LRConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
		virtual void LRElasticConstitutiveLaw(MPMBase *,Matrix3 &,Matrix3 &,Matrix3 &,Matrix3 &,Matrix3 &,int,void *,ResidualStrains *) const;
		virtual void SRConstitutiveLaw2D(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
		virtual void SRConstitutiveLaw3D(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
		virtual ElasticProperties *GetElasticPropertiesPointer(void *) const;

		// methods for softening materials
		virtual double GetAcOverVp(int,MPMBase *,Vector *) const;
		virtual bool SoftenAxis(double,double *,int,int,double,double,SofteningLaw *,
								  double,double,double,double &,double *,double &,bool &) const;
		virtual bool SoftenTwoAxes(double *,double,double,double,double,double,SofteningLaw *,int,int,
			double,double,double,double,double,SofteningLaw *,int,int,double &,double &,double &,bool &) const;
		virtual void PostFailureUpdate(double &,double &,double &,Tensor *,Tensor *,Tensor *,Matrix3,
									   double,double,double,double,bool) const;
		virtual void UpdateCrackingStrainStress(int,Tensor *,Tensor *,double,double,double,Tensor *,Tensor *,Matrix3) const;
#else
        virtual double GetStressStrainZZ(double,double,double,double,double,int);
#endif

	protected:
		double prop1,prop2;
#ifdef MPM_CODE
		ElasticProperties pr;
		double Cadota;
		int useLargeRotation;				// Hypoelastic material with large rotation methods
#else
        double prop3;
#endif

		// methods
        virtual const char *SetAnalysisProps(int,double,double,double,double,double,
                        double,double,double,double,double,double,double,double,double,double);
};

#endif
