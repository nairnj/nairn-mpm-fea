/********************************************************************************
    OrthoPlasticSoftening.hpp
    nairn-mpm-fea
    
    Created by John Nairn, Oct 3, 2020.
    Copyright (c) 2020 John A. Nairn, All rights reserved.

	Dependencies
        OrthoSofening.hpp TransIsoSoftening.hpp TransIsotropic.hpp
        Elastic.hpp MaterialBase.hpp
********************************************************************************/

#ifndef ORTHOPLASTICSOFTENING

#define ORTHOPLASTICSOFTENING 56

#include "Materials/OrthoSoftening.hpp"
#include "Materials/AnisoPlasticity.hpp"

class OrthoPlasticSoftening : public OrthoSoftening
{
    public:
        // constructors and destructors
		OrthoPlasticSoftening(char *,int);
		
		// initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;
        virtual char *InitHistoryData(char *,MPMBase *);
		virtual void ResetHistoryData(char *,MPMBase *);
        virtual int NumberOfHistoryDoubles(void) const;

		// contitutive law methods
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int,Tensor *) const;
        virtual void UpdateCrackingStrain(int,Tensor *,double,double,double,Matrix3,double *) const;
        virtual void LoadCrackHillProperties(int,HillProperties *,CrackAxisProperties *d,int) const;
        virtual double GetYield(AnisoPlasticProperties *) const;
        virtual double GetGPrime(AnisoPlasticProperties *) const;
        virtual Tensor SolveForPlasticIncrement(MPMBase *,int,double,Tensor &,AnisoPlasticProperties *,
                                                CrackAxisProperties *d,HillProperties *h) const;
        virtual void GetDfCdf(Tensor &,int,AnisoPlasticProperties *,CrackAxisProperties *,HillProperties *) const;
        virtual ElasticProperties *GetElasticPropertiesPointer(void *) const;

		// accessors
		virtual const char *MaterialType(void) const;
        virtual int SizeOfMechanicalProperties(int &) const;
        virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *,int) const;
        virtual double *GetSoftHistoryPtr(MPMBase *mptr) const;
        virtual int AltStrainContains(void) const;
        virtual bool GetCrackingStrain(MPMBase *,Tensor *,bool,Matrix3 *) const;
        virtual Vector GetCrackingCOD(MPMBase *,bool) const;

    protected:
        // Harding constantsL Y/Yref = 1 + Khard alpha^nhard
        double Khard,nhard,exphard,Kexp,alphaMax,hmax,slopeMin;
        int hardStyle;
        double syxx,syyy,syzz,tyyz,tyxz,tyxy;
        double sqrt23OversigmaYref;
        HillProperties hMAS;        // in material axis system

};

#endif

