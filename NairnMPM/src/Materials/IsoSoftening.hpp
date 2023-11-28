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

// vector 3D
typedef struct {
    double deltaN;
    double deltaS;
    double damageN;
    double damageS;
    double decxx;
    double dgcxy;
    double dgcxz;
    double ecxx;
} DamageState;

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
		virtual bool AcceptInitiationLaw(FailureSurface *,int);
		virtual bool AcceptSofteningLaw(SofteningLaw *,int);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;
	
		// history data
		virtual char *InitHistoryData(char *,MPMBase *);
		virtual void ResetHistoryData(char *,MPMBase *);
   		virtual int NumberOfHistoryDoubles(void) const;
        virtual void SetInitialConditions(InitialCondition *,int,bool);
		virtual Vector GetDamageNormal(MPMBase *,bool) const;
	
		// methods
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int,Tensor *) const;
		virtual void DamageEvolution(MPMBase *,int,double *,Tensor &,Tensor &,double,double,ResidualStrains *,
								  			Matrix3 &,Matrix3 &,ElasticProperties *,Tensor *,double,double) const;
		virtual bool GetRToCrack(Matrix3 *,double *,bool,int) const;
	
		// isotropic elasticity methods
		virtual bool CoupledShearSoftening(double,double,double *,Tensor &,double,double,double,double,
											double,double &,double &,double &,bool &) const;
		virtual bool OvoidSoftening(MPMBase *,bool,double,double,double,DamageState *,
                                    Tensor &,double,double,double,double,double,
									double,double,double,double &,double &,bool &) const;
		virtual bool CoupledCuboidSoftening(MPMBase *,bool,double,double,double,DamageState *,
									Tensor &,double,double,double,double,double,
											double,double,double,double &,double &,bool &) const;
		virtual Tensor GetStressIncrement(Tensor &,int,void *) const;
		virtual void AcceptTrialStress(MPMBase *,Tensor &,Tensor *,int,Matrix3 *,void *,Tensor &,double,double) const;
	
		// accessors
		virtual const char *MaterialType() const;
		virtual int AltStrainContains(void) const;
		virtual bool GetCrackingStrain(MPMBase *,Tensor *,bool,Matrix3 *) const;
		virtual Vector GetCrackingCOD(MPMBase *,bool) const;
		virtual double *GetSoftHistoryPtr(MPMBase *) const;
		virtual int GetTractionFailureSurface(void) const;
        virtual void SetRelativeStrength(MPMBase *,double);
        virtual void SetRelativeToughness(MPMBase *,double);

	protected:
		FailureSurface *initiationLaw;
		SofteningLaw *softeningModeI;
		SofteningLaw *softeningModeII;
		double en0,gs0;
		int tractionFailureSurface;
		double softenCV,wAlpha,wV0,wGam1A;
        int softenStatsMode;
		int distributionMode;
		double frictionCoeff;
        double pdOvoidTolerance;
        int maxOvoidPasses;
};

#endif
