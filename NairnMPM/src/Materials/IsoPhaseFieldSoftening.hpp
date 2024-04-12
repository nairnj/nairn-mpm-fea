/********************************************************************************
	IsoPhaseFieldSoftening.hpp
	nairn-mpm-fea
 
	Created by John Nairn, Sept. 27, 2021
	Copyright (c) 2021 John Nairn, All rights reserved.
 
	Dependencies
		IsotropicMat.hpp (Elastic.hpp MaterialBase.hpp)
 ********************************************************************************/

#ifndef ISOPHASEFIELDSOFTENING

#define ISOPHASEFIELDSOFTENING 57

#include "Materials/IsotropicMat.hpp"

class PhaseFieldDiffusion;

enum { GEN_QUADRATIC=0, FINITE_EXPONENTIAL, PF_LINEAR_SOFTENING };

enum { TENSILE_EIGENVALUES=0, PRESSURE_SHEAR };

enum { HISTORY_VALUE=0, PHASE_FAILURE_STATE, PHASE_SOLVED, DELTA_PHASE_SOLVED };

class IsoPhaseFieldSoftening : public IsotropicMat
{
	public:

		// constructors and destructors
		IsoPhaseFieldSoftening(char *,int);
	
		// initialize
		virtual char *InputMaterialProperty(char *,int &input,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;

		// history data
		virtual char *InitHistoryData(char *pchr,MPMBase *mptr);
		virtual int NumberOfHistoryDoubles(void) const;
		virtual void SetInitialConditions(InitialCondition *,int,bool);
	
		// methods
		virtual void LRConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
	
		// accessors
        virtual const char *MaterialType() const;
        virtual double GetDegradation(MPMBase *,double &) const;

        // phase field diffusion accessors
        virtual bool SupportsPhaseField(int) const;
        virtual double MaximumPhaseDiffusivity(int) const;
        virtual Tensor GetPhaseFieldDTensor(int,MPMBase *) const;
        virtual double GetPhaseFieldViscosity(int) const;
		virtual bool GetParticleDiffusionSource(DiffusionTask *,MPMBase *,double,double *) const;
        virtual void AdjustPhaseFieldValue(DiffusionTask *,double,double,double &,double &,double &,double) const;
        virtual void StorePhaseFieldValue(int,MPMBase *,double) const;
        virtual void StorePhaseFieldDelta(int,MPMBase *,double) const;

	protected:
		double viscosity;
		double ell;
		int gdMode;
		double kStability;
		bool hasStability;
		double garg;
		double ec1,ec2,ec3;
		int partition;
		int gargPicked;
		double Psii;
        PhaseFieldDiffusion *phaseTask;
		int taskNum;
		Tensor pfDTensor;
};

#endif
