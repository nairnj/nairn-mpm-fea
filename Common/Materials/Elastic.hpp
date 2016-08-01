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
// alpha and beta are thermal and moisture expansion coefficients
// although some elements are used for other things
typedef struct {
	double C[6][6];
	double alpha[8];
	double beta[6];
} ElasticProperties;
#endif

class Elastic : public MaterialBase
{
    public:
        
        // constructors and destructors
        Elastic();
        Elastic(char *);
		
		// initialize
		char *InputMaterialProperty(char *,int &,double &);
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
		virtual void LRElasticConstitutiveLaw(MPMBase *,Matrix3,Matrix3,Matrix3,Matrix3,Matrix3 *,int,void *,ResidualStrains *) const;
		virtual void SRConstitutiveLaw2D(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
		virtual void SRConstitutiveLaw3D(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
		virtual ElasticProperties *GetElasticPropertiesPointer(void *) const;
		virtual double GetAcOverVp(int,MPMBase *,Vector *) const;
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
