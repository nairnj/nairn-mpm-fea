/********************************************************************************
    IsoPlasticity.hpp
    nairn-mpm-fea
    
    Created by John Nairn, June 16, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		IsotropicMat.hpp (Elastic.hpp, MaterialBase.hpp)
********************************************************************************/

#ifndef _ISOPLASTICITY_

#define _ISOPLASTICITY_

#define ISOPLASTICITY 9

// plastic law properties
typedef struct {
	double Gred;
	double Kred;
	double psRed;
	double psLr2G;
	double psKred;
	void *hardProps;
	double delVLowStrain;
	double QAVred;
} PlasticProperties;

class HardeningLawBase;

#include "Materials/IsotropicMat.hpp"

class IsoPlasticity : public IsotropicMat
{
    public:
        // constructors and destructors
		IsoPlasticity();
		IsoPlasticity(char *matName);
		
		// initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
        virtual bool AcceptHardeningLaw(HardeningLawBase *,int );
	
		// history data
		virtual char *InitHistoryData(char *,MPMBase *);
		virtual double GetHistory(int,char *) const;
 	
		// const methods
        virtual void PrintMechanicalProperties(void) const;
		
		// methods
        virtual int SizeOfMechanicalProperties(int &) const;
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *,int) const;
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int) const;
		virtual void LRConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
        virtual void LRPlasticityConstLaw(MPMBase *,double,double,double,double,double,int,
                                        double,double,PlasticProperties *,ResidualStrains *,Matrix3 *) const;
        virtual void LRPlasticityConstLaw(MPMBase *,double,double,double,double,double,
                                        double,double,int,double,double,PlasticProperties *,ResidualStrains *,Matrix3 *) const;
		virtual void SRConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
		virtual void SRPlasticityConstLaw2D(MPMBase *,Matrix3,double,int,
										double,double,PlasticProperties *,ResidualStrains *,Matrix3 *) const;
		virtual void SRPlasticityConstLaw3D(MPMBase *,Matrix3,double,int,
										double,double,PlasticProperties *,ResidualStrains *,Matrix3 *) const;
		
		// custom methods: Find yield function and solve for lambda
		virtual void UpdatePressure(MPMBase *,double,int,PlasticProperties *,ResidualStrains *,double,double &,double &) const;
        virtual double GetMagnitudeSFromDev(Tensor *,int) const;
		virtual void GetDfDsigma(double,Tensor *,int,Tensor *) const;
		
		// accessors
        virtual Tensor GetStress(Tensor *,double,MPMBase *) const;
       int MaterialTag(void) const;
        const char *MaterialType(void) const;
		virtual int AltStrainContains(void) const;
		
    protected:
		PlasticProperties pr;
		double G0red;
        HardeningLawBase *plasticLaw;

};

#endif

