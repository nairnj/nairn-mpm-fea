/********************************************************************************
    TransportTask.hpp
    nairn-mpm-fea
    
    Created by John Nairn on 7/18/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _TRANSPORTTASK_

#define _TRANSPORTTASK_

class NodalPoint;
class MPMBase;

class TransportTask
{
    public:
        TransportTask *nextTask;
       
        // constructors and destructors
        TransportTask();
        virtual ~TransportTask();
        
        // overridable methods
		virtual TransportTask *TransportOutput(void);
		
		// pure virtual required methods
		
		// Return name of this task
		virtual const char *TaskName(void) = 0;
		
		// adjust time for given cell size if needed
		virtual TransportTask *TransportTimeStep(int,double,double *) = 0;
		
		// Task 1 Extrapolation of transport property to the grid
		virtual TransportTask *Task1Extrapolation(NodalPoint *,MPMBase *,double) = 0;
		virtual TransportTask *Task1Reduction(NodalPoint *,NodalPoint *) = 0;

		// Get grid values
		virtual TransportTask *GetNodalValue(NodalPoint *) = 0;
		
		// impose grid-based BCs
		virtual void ImposeValueBCs(double) = 0;
	
		// Get transport property gradient on the particles
		virtual TransportTask *GetGradients(double) = 0;
		
		// find forces for transport calculation
		virtual TransportTask *AddForces(NodalPoint *,MPMBase *,double,double,double,double,TransportProperties *) = 0;
        virtual TransportTask *CopyForces(NodalPoint *,NodalPoint *) = 0;

		// adjust forces at grid points with transport BCs
		virtual TransportTask *SetTransportForceBCs(double) = 0;
		
		// find transport rates on the nodes
		virtual TransportTask *TransportRates(NodalPoint *,double) = 0;
		
		// increment transport rate
		virtual TransportTask *IncrementTransportRate(const NodalPoint *,double,double &) const = 0;
		
		// increment particle transpoprt value
		virtual TransportTask *MoveTransportValue(MPMBase *,double,double) const = 0;
		
		// if needed for SZS or USAVG, update value on the grid
		virtual TransportTask *UpdateNodalValues(double) = 0;
		
		// increment transport rate when updating particle strain
		virtual double IncrementValueExtrap(NodalPoint *,double) const = 0;
		
		// find change in transport value on the particle from grid results
		virtual double GetDeltaValue(MPMBase *,double) const = 0;

};

// globals
extern TransportTask *transportTasks;

#endif

