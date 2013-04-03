/********************************************************************************
    IsoPlasticity.hpp
    NairnMPM
    
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
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int);
		virtual void DeleteCopyOfMechanicalProps(void *,int) const;
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *);
        virtual void PlasticityConstLaw(MPMBase *,double,double,double,double,double,double,int,
                                        double,double,double,PlasticProperties *,ResidualStrains *);
        virtual void PlasticityConstLaw(MPMBase *,double,double,double,double,double,double,double,double,
                                        double,double,int,double,double,double,PlasticProperties *,ResidualStrains *);
		
		// custom methods: Find yield function and solve for lambda
		virtual void UpdatePressure(MPMBase *,double &,double,int,PlasticProperties *,ResidualStrains *);
        virtual double GetMagnitudeSFromDev(Tensor *,int);
		virtual void GetDfDsigma(double,Tensor *,int);
		virtual void ElasticUpdateFinished(MPMBase *,int,double);
		
		// accessors
        virtual Tensor GetStress(Tensor *,double) const;
		virtual double GetHistory(int,char *) const;
        virtual bool PartitionsElasticAndPlasticStrain(void);
        int MaterialTag(void) const;
        const char *MaterialType(void) const;
		
    protected:
		PlasticProperties pr;
		double dfdsxx,dfdsyy,dfdtxy,dfdszz,dfdtxz,dfdtyz;
        double delVLowStrain;
		double G0red;
    
        HardeningLawBase *plasticLaw;

};

#endif

