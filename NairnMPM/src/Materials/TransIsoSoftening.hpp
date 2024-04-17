/********************************************************************************
	TransIsoSoftening.hpp
	nairn-mpm-fea

	Created by John Nairn on 4 AUG 2016
	Copyright (c) 2015 John A. Nairn, All rights reserved.

	Dependencies
		TransIsotropic.hpp (Elastic.hpp, MaterialBase.hpp)
********************************************************************************/

#ifndef TRANSISOSOFTENING1

#define TRANSISOSOFTENING1 51
#define TRANSISOSOFTENING2 52

#include "Materials/TransIsotropic.hpp"

typedef struct {
	double C12C11;
	double C13C11;
	double C11;
	double C12;
	double C13;
	double C22;
	double C23;
	double C33;
	double C44;
	double C55;
	double C66;
	SofteningLaw *fnLaw;
	SofteningLaw *fXYLaw;
	SofteningLaw *fXZLaw;
	double sigNc;
	double tauXYc;
	double tauXZc;
	double vzxc;
	double vzyc;
} CrackAxisProperties;

class FailureSurface;
class SofteningLaw;
class InitialCondition;

// Damage processes
enum { DZ_DAMAGE=0,DX_DAMAGE,DY_DAMAGE };

extern double predamageState,damageState;

class TransIsoSoftening : public TransIsotropic
{
	public:
	
		// constructors and destructors
		TransIsoSoftening(char *,int);
		
		// initialize
		virtual char *InputMaterialProperty(char *,int &,double &);
		virtual bool AcceptInitiationLaw(FailureSurface *,int);
		virtual bool AcceptSofteningLaw(SofteningLaw *,int);
        void SwapLaws(SofteningLaw **,SofteningLaw **);
 		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;
	
		// history data
		virtual char *InitHistoryData(char *,MPMBase *);
		virtual void ResetHistoryData(char *,MPMBase *);
		virtual void VerifyStability(double,double,double,SofteningLaw *,const char *,const char *) const;
  		virtual int NumberOfHistoryDoubles(void) const;
        virtual void SetInitialConditions(InitialCondition *,int,bool);
		virtual Vector GetDamageNormal(MPMBase *,bool) const;
        virtual int GetDForm(double) const;

		// methods
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int,Tensor *) const;
        void DamageEvolution(MPMBase *,int,double *,Tensor &,Tensor &,Tensor &,
                             CrackAxisProperties &,Matrix3 &,Matrix3 &,Tensor &,int,double,double) const;
		virtual int DecodeDamageInitiation(int,Vector *,int,double *) const;
		virtual void LoadCrackAxisProperties(int,CrackAxisProperties *,int,ElasticProperties *) const;
		virtual bool GetRToCrack(Matrix3 *,double *,bool,int) const;
        virtual Tensor GetAnisoResStrains(double &,double &,double &,ElasticProperties *,ResidualStrains *,int) const;
        bool OvoidSoftening(MPMBase *mptr,bool is2D,double *soft,
                            double den,double dgxy,double dgxz,Tensor &str,CrackAxisProperties *d,
                            double sigmaN,double scaleN,double sigmaXY,double scaleXY,double sigmaXZ,double scaleXZ,
                            double &dispNEnergy,double &dispXYEnergy,double &dispXZEnergy,bool &criticalStrain,
                            double,double &,double &,double &) const;

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
		SofteningLaw *softeningAI;          // XX for ortho
		SofteningLaw *softeningI;           // YY for ortho
		SofteningLaw *softeningTII;         // XYY for ortho
		SofteningLaw *softeningAII;         // XYX for ortho
		SofteningLaw *softeningII;          // YZY for ortho
		int tractionFailureSurface;
		double softenCV,wAlpha,wV0,wGam1A;
		int softenStatsMode;
		int distributionMode;
 		double frictionCoeff;
};

#endif
