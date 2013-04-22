/********************************************************************************
    ConductionTask.hpp
    NairnMPM
    
    Created by John Nairn on Fri Oct 15 2004
    Copyright (c) 2004 John A. Nairn, All rights reserved.
	
	Dependencies
		TransportTask.hpp
********************************************************************************/

#ifndef _CONDUCTIONTASK_

#define _CONDUCTIONTASK_

#include "Custom_Tasks/TransportTask.hpp"

class MPMBase;
class NodalPoint;
class CrackSegment;

class ConductionTask : public TransportTask
{
    public:
		static bool active,crackTipHeating,energyCoupling,AVHeating;
		
        // constructors and destructors
		ConductionTask();
        
        // standard methods
		virtual const char *TaskName(void);
		virtual TransportTask *TransportOutput(void);
		virtual TransportTask *TransportTimeStep(int,double,double *);
		virtual TransportTask *Task1Extrapolation(NodalPoint *,MPMBase *,double);
		virtual TransportTask *Task1Reduction(NodalPoint *,NodalPoint *);
		virtual TransportTask *GetNodalValue(NodalPoint *);
		virtual void ImposeValueBCs(double);
		virtual TransportTask *GetGradients(double);
		virtual TransportTask *AddForces(NodalPoint *,MPMBase *,double,double,double,double,TransportProperties *);
		virtual TransportTask *CopyForces(NodalPoint *,NodalPoint *);
		virtual TransportTask *SetTransportForceBCs(double);
		virtual TransportTask *TransportRates(NodalPoint *,double);
		virtual TransportTask *IncrementTransportRate(NodalPoint *,double,double &) const;
		virtual TransportTask *MoveTransportValue(MPMBase *,double,double) const;
		virtual TransportTask *UpdateNodalValues(double);
		virtual double IncrementValueExtrap(NodalPoint *,double) const;
		virtual double GetDeltaValue(MPMBase *,double) const;
		
		// custom methods
		void AddCrackTipHeating(void);
		void StartCrackTipHeating(CrackSegment *,Vector &,double);
    
        // class methods
        static void ThermodynamicsOutput(void);
        static bool IsSystemIsolated(void);
		
    private:
};

// globals
extern ConductionTask *conduction;

#endif
