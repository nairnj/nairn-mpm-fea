/********************************************************************************
    DiffusionTask.hpp
    NairnMPM
    
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

class DiffusionTask : public TransportTask
{
    public:
		static bool active;
		static double reference;
		
        // constructors and destructors
		DiffusionTask();
        
        // standard methods
		virtual TransportTask *TransportOutput(void);
		virtual const char *TaskName(void);
		virtual TransportTask *TransportTimeStep(int,double,double *);
		virtual TransportTask *Task1Extrapolation(NodalPoint *,MPMBase *,double);
		virtual TransportTask *GetNodalValue(NodalPoint *);
		virtual void ImposeValueBCs(double);
		virtual TransportTask *GetGradients(double);
		virtual TransportTask *AddForces(double *,MPMBase *,double,double,double,double,TransportProperties *);
		virtual TransportTask *AddForcesFromBuffer(NodalPoint *,double);
		virtual TransportTask *SetTransportForceBCs(double);
		virtual TransportTask *TransportRates(double);
		virtual TransportTask *IncrementTransportRate(NodalPoint *,double,double &) const;
		virtual TransportTask *MoveTransportValue(MPMBase *,double,double) const;
		virtual TransportTask *UpdateNodalValues(double);
		virtual double IncrementValueExtrap(NodalPoint *,double) const;
		virtual double GetDeltaValue(MPMBase *,double) const;
		
    private:
};

// globals
extern DiffusionTask *diffusion;

#endif
