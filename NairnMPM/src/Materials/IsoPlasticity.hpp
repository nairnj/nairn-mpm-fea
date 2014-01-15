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
        virtual char *InputMat(char *,int &);
        virtual const char *VerifyAndLoadProperties(int);
        virtual void SetHardeningLaw(char *);
		virtual char *InitHistoryData(void);
	
		// const methods
        virtual void PrintMechanicalProperties(void) const;
		
		// methods
        virtual int SizeOfMechanicalProperties(int &) const;
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *) const;
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
        virtual void PlasticityConstLaw(MPMBase *,double,double,double,double,double,double,int,
                                        double,double,double,PlasticProperties *,ResidualStrains *) const;
        virtual void PlasticityConstLaw(MPMBase *,double,double,double,double,double,double,double,double,
                                        double,double,int,double,double,double,PlasticProperties *,ResidualStrains *) const;
		
		// custom methods: Find yield function and solve for lambda
		virtual void UpdatePressure(MPMBase *,double &,double,int,PlasticProperties *,ResidualStrains *,double) const;
        virtual double GetMagnitudeSFromDev(Tensor *,int) const;
		virtual void GetDfDsigma(double,Tensor *,int,Tensor *) const;
		virtual void ElasticUpdateFinished(MPMBase *,int,double) const;
		
		// accessors
        virtual Tensor GetStress(Tensor *,double) const;
		virtual double GetHistory(int,char *) const;
        virtual bool PartitionsElasticAndPlasticStrain(void);
        int MaterialTag(void) const;
        const char *MaterialType(void) const;
		
    protected:
		PlasticProperties pr;
		double G0red;
        HardeningLawBase *plasticLaw;

};

#endif

