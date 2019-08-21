/********************************************************************************
	IsoSoftening.hpp
	nairn-mpm-fea

	Created by John Nairn on Jan 26, 2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.

	Dependencies
		IsotropicMat.hpp (Elastic.hpp MaterialBase.hpp)
********************************************************************************/

#ifndef ISOSOFTENING

#define ISOSOFTENING 50

#include "Materials/IsotropicMat.hpp"

class FailureSurface;
class SofteningLaw;
class InitialCondition;

extern double predamageState,damageState;

class IsoSoftening : public IsotropicMat
{
	public:
	
		// constructors and destructors
		IsoSoftening(char *,int);
		
		// initialize
		virtual char *InputMaterialProperty(char *,int &,double &);
		bool AcceptInitiationLaw(FailureSurface *,int);
		bool AcceptSofteningLaw(SofteningLaw *,int,int);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;
	
		// history data
		virtual char *InitHistoryData(char *,MPMBase *);
   		virtual int NumberOfHistoryDoubles(void) const;
        virtual void SetInitialConditions(InitialCondition *,MPMBase *,bool);
		virtual Vector GetDamageNormal(MPMBase *,bool) const;
	
		// methods
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int,Tensor *) const;
		virtual Vector DamageEvolution(MPMBase *,int,double *,Tensor &,Tensor &,double,double,ResidualStrains *,
								  			Matrix3 &,Matrix3 &,ElasticProperties *,Tensor *,double) const;
		virtual bool GetRToCrack(Matrix3 *,double *,int,int) const;
	
		// isotropic elasticity methods
		virtual Tensor GetStressIncrement(Tensor &,int,void *) const;
		virtual void AcceptTrialStress(MPMBase *,Tensor &,Tensor *,int,Matrix3 *,void *,Tensor &,double,double) const;
	
		// accessors
		virtual const char *MaterialType() const;
		virtual int AltStrainContains(void) const;
		virtual double *GetSoftHistoryPtr(MPMBase *) const;
	
	protected:
		FailureSurface *initiationLaw;
		SofteningLaw *softeningModeI;
		SofteningLaw *softeningModeII;
		double en0,gs0;
		int shearFailureSurface;
		double softenCV;
        int softenCVMode;
};

#endif
