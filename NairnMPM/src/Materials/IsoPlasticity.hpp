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
        virtual const char *VerifyProperties(int);
        virtual void SetHardeningLaw(char *);
		virtual void InitialLoadMechProps(int,int);
        virtual void PrintMechanicalProperties(void);
		virtual char *InitHistoryData(void);
		
		// methods
        virtual void LoadMechanicalProps(MPMBase *,int);
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int);
        virtual void PlasticityConstLaw(MPMBase *,double,double,double,double,double,double,int,
                                        double,double,double);
        virtual void PlasticityConstLaw(MPMBase *,double,double,double,double,double,double,double,double,
                                        double,double,int,double,double,double);
		
		// custom methods: Find yield function and solve for lambda
		virtual void UpdatePressure(MPMBase *,double &,double,int);
        virtual double GetMagnitudeSFromDev(Tensor *,int);
		virtual void GetDfDsigma(double,Tensor *,int);
		virtual void ElasticUpdateFinished(MPMBase *,int,double);
		
		// accessors
        virtual Tensor GetStress(Tensor *,double);
		virtual double GetHistory(int,char *);
        virtual bool HasPlasticStrainForGradient(void);
        int MaterialTag(void);
        const char *MaterialType(void);
		
    protected:
		double Gred,Kred;
		double psRed,psLr2G,psKred;
		double dfdsxx,dfdsyy,dfdtxy,dfdszz,dfdtxz,dfdtyz;
        double delVLowStrain;
    
        HardeningLawBase *plasticLaw;

};

#endif

