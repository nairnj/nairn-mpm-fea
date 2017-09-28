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

class MPMBase;
class NodalPoint;
class MatVelocityField;

class DiffusionTask : public TransportTask
{
    public:
		static bool active;
		static double reference;

		// initialize
		virtual TransportTask *Initialize(void);
	
		// mass and momentum
		virtual void ZeroTransportGradients(MPMBase *);
		virtual void AddTransportGradients(MPMBase *,Vector *,NodalPoint *,short);
	
		// grid forces
		virtual TransportTask *AddForces(NodalPoint *,MPMBase *,double,double,double,double,TransportProperties *,short,int);
	
		// update particles task
		virtual TransportTask *MoveTransportValue(MPMBase *,double,double) const;
	
        // accessors
        virtual const char *TaskName(void);
        virtual TransportTask *TransportTimeStepFactor(int,double *);
        virtual double GetTransportMassAndValue(MPMBase *,double *);
        virtual TransportField *GetTransportFieldPtr(NodalPoint *) const;
        virtual NodalValueBC *GetFirstBCPtr(void) const;
		virtual MatPtLoadBC *GetFirstFluxBCPtr(void) const;
        virtual double *GetParticleValuePtr(MPMBase *mptr) const;
        virtual double *GetPrevParticleValuePtr(MPMBase *mptr) const;
	
    private:
};

// globals
extern DiffusionTask *diffusion;

#endif
