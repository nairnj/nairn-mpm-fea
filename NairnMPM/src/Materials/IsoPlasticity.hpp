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
		IsoPlasticity(char *,int);
		
		// initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
        virtual bool AcceptHardeningLaw(HardeningLawBase *,int );
		virtual HardeningLawBase *GetPlasticLaw(void) const;

		// history data
		virtual char *InitHistoryData(char *,MPMBase *);
   		virtual int NumberOfHistoryDoubles(void) const;
 	
		// const methods
        virtual void PrintMechanicalProperties(void) const;
		
		// methods
        virtual int SizeOfMechanicalProperties(int &) const;
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *,int) const;
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int) const;
        virtual void PlasticityConstLaw(MPMBase *,Matrix3,double,int,double,double,
                                        PlasticProperties *,ResidualStrains *,Matrix3 *,bool) const;
		
		// custom methods: Find yield function and solve for lambda
		virtual void UpdatePressure(MPMBase *,double,int,PlasticProperties *,ResidualStrains *,double,double,double &,double &) const;
		virtual void GetDfDsigma(double,Tensor *,int,Tensor *) const;
		virtual double GetPlasticPotential(Tensor *,MPMBase *,int,double,HardeningAlpha *,void *) const;
		virtual double RRPlasticIncrement(Tensor *,MPMBase *,int,double,double,HardeningAlpha *,PlasticProperties *) const;
	
		// accessors
        virtual Tensor GetStress(Tensor *,double,MPMBase *) const;
		virtual void SetStress(Tensor *,MPMBase *) const;
		virtual void IncrementThicknessStress(double,MPMBase *) const;
        const char *MaterialType(void) const;
		virtual int AltStrainContains(void) const;
		virtual bool SupportsArtificialViscosity(void) const;

    protected:
		PlasticProperties pr;
		double G0red;
        HardeningLawBase *plasticLaw;

};

#endif

