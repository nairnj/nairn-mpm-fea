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

enum { NO_DIFFUSION=0,MOISTURE_DIFFUSION,POROELASTICITY_DIFFUSION };

class MPMBase;
class NodalPoint;
class MatVelocityField;

class DiffusionTask : public TransportTask
{
    public:
		static int active;
		static double reference;
		static double viscosity;

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
		virtual TransportTask *MoveTransportValue(MPMBase *,double,double) const;
	
        // accessors
        virtual const char *TaskName(void);
        virtual TransportTask *TransportTimeStepFactor(int,double *);
        virtual TransportField *GetTransportFieldPtr(NodalPoint *) const;
        virtual NodalValueBC *GetFirstBCPtr(void) const;
		virtual MatPtLoadBC *GetFirstFluxBCPtr(void) const;
        virtual double *GetParticleValuePtr(MPMBase *mptr) const;
        virtual double *GetPrevParticleValuePtr(MPMBase *mptr) const;
	
		// static methods
		static bool HasDiffusion(void);
		static bool HasPoroelasticity(void);
		static bool HasFluidTransport(void);
		static double RescalePotential(void);
		static double RescaleFlux(void);
	
    private:
};

// globals
extern DiffusionTask *diffusion;

#endif
