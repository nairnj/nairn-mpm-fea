/********************************************************************************
    DiffusionTask.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Mon Mar 08 2004
    Copyright (c) 2003 John A. Nairn, All rights reserved.
	
	Dependencies
		TransportTask.hpp
********************************************************************************/

#ifndef _DIFFUSIONTASK_

#define _DIFFUSIONTASK_

#include "Custom_Tasks/TransportTask.hpp"

// Move to MaterialBase to be more accessible
//enum { NO_DIFFUSION=0,MOISTURE_DIFFUSION,POROELASTICITY_DIFFUSION,
//    FRACTURE_PHASE_FIELD,BATTERY_PHASE_FIELD,CONDUCTION_PHASE_FIELD,NUM_DUFFUSION_OPTIONS };

class MPMBase;
class NodalPoint;
class MatVelocityField;

class DiffusionTask : public TransportTask
{
    public:
		double reference;
		double viscosity;
		int style;
        bool noLimit;

		// constructors and destructors
		DiffusionTask(double,double,int);

		// initialize
		virtual TransportTask *Initialize(void);
	
		// mass and momentum
		virtual TransportTask *Task1Extrapolation(NodalPoint *,MPMBase *,double,short,int);
		virtual double GetVpCTp(MPMBase *);
		virtual void ZeroTransportGradients(MPMBase *);
		virtual void AddTransportGradients(MPMBase *,Vector *,NodalPoint *,short);

		// grid forces
		virtual TransportTask *AddForces(NodalPoint *,MPMBase *,double,double,double,double,TransportProperties *,short,int);
	
		// update particles task
		virtual void AdjustRateAndValue(MPMBase *,double &,double &,double &,double) const;
		virtual void GetDeltaValue(MPMBase *,double,double *) const;
		virtual TransportTask *MoveTransportValue(MPMBase *,double,double,double) const;
	
        // accessors
        virtual const char *TaskName(void);
		virtual const char *StyleName(void);
        virtual TransportTask *TransportTimeStepFactor(int,double *);
        virtual TransportField *GetTransportFieldPtr(NodalPoint *) const;
        virtual NodalValueBC *GetFirstBCPtr(void) const;
		virtual MatPtLoadBC *GetFirstFluxBCPtr(void) const;
		virtual double GetParticleValue(MPMBase *) const;
		virtual double *GetParticleValuePtr(MPMBase *) const;
		virtual double GetPrevParticleValue(MPMBase *) const;
        virtual double *GetPrevParticleValuePtr(MPMBase *) const;
		virtual double GetDeltaConcentration(MPMBase *) const;
		virtual int GetNumber(void) const;
		virtual void SetNumber(int);
		virtual int GetStyle(void) const;
	
		// static methods
		static bool HasDiffusion(void);
#ifdef POROELASTICITY
		static bool HasPoroelasticity(void);
#endif
		static bool HasFluidTransport(void);
		static double RescalePotential(int);
		static double RescaleFlux(void);
		static DiffusionTask *FindDiffusionTask(int);
		static int FindDiffusionNumber(int);
		static void CountDiffusionTasks(void);
		static double GetMinimumTimeStep(void);
		static void SetDiffusionXPIC(bool);

    protected:
		int number;
};

// globals
extern DiffusionTask *diffusion;
extern DiffusionTask *otherDiffusion;
extern int numDiffusion;

#endif
