/********************************************************************************
	IsoPlasticSoftening.hpp
	nairn-mpm-fea
 
	Created by John Nairn on Oct 9, 2017
	Copyright (c) 2017 John A. Nairn, All rights reserved.
 
	Dependencies
		IsoSoftening.hpp IsotropicMat.hpp (Elastic.hpp MaterialBase.hpp)
 ********************************************************************************/

#ifndef ISOPLASTICSOFTENING

#define ISOPLASTICSOFTENING 53

#include "Materials/IsoSoftening.hpp"

// plastic law properties
typedef struct {
	void *hardProps;
	void *elasticProps;
} SofteningPlasticProperties;


class IsoPlasticSoftening : public IsoSoftening
{
	public:
	
		// constructors and destructors
		IsoPlasticSoftening(char *,int);
	
		// initialize
		virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual bool AcceptHardeningLaw(HardeningLawBase *,int );
		virtual HardeningLawBase *GetPlasticLaw(void) const;

		// const methods
		virtual void PrintMechanicalProperties(void) const;
	
    	// history data
    	virtual char *InitHistoryData(char *,MPMBase *);
   		virtual int NumberOfHistoryDoubles(void) const;
    
		// methods
		virtual int SizeOfMechanicalProperties(int &) const;
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *,int) const;
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int,Tensor *) const;
		virtual void UpdateCrackingStrain(int,Tensor *,double,double,double,Matrix3,double *) const;
	
		// accessors
		const char *MaterialType(void) const;
		virtual int AltStrainContains(void) const;
		virtual bool GetCrackingStrain(MPMBase *,Tensor *,bool,Matrix3 *) const;
		virtual Vector GetCrackingCOD(MPMBase *,bool) const;
		virtual double *GetSoftHistoryPtr(MPMBase *) const;
#ifdef SMOOTHED_STRESS
		virtual double GetNonLocalRadius(void) const;
#endif
	
	protected:
		HardeningLawBase *plasticLaw;
    	int softHistoryOffset;
		double Kred,Gred,kappa;
};

#endif
