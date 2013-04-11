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
		static bool active,crackTipHeating,energyCoupling, frictionalHeating, hardCodedHeat; //modiftf #frictionalheating #hardcodedheat
		static double dTemperature;
		
        // constructors and destructors
		ConductionTask();
        
        // standard methods
		virtual const char *TaskName(void);
		virtual TransportTask *TransportOutput(void);
		virtual TransportTask *TransportTimeStep(int,double,double *);
		virtual TransportTask *Task1Extrapolation(NodalPoint *,MPMBase *,double);
		virtual void GetValues(double);
		virtual void GetGradients(double);
		virtual TransportTask *AddForces(NodalPoint *,MPMBase *,double,double,double,double);
		virtual TransportTask *SetTransportForceBCs(double);
		virtual TransportTask *TransportRates(double);
		virtual TransportTask *IncrementTransportRate(NodalPoint *,double);
		virtual TransportTask *MoveTransportValue(MPMBase *,double);
		virtual TransportTask *UpdateNodalValues(double);
		virtual TransportTask *IncrementValueExtrap(NodalPoint *,double);
		virtual TransportTask *GetDeltaValue(MPMBase *);
		
		// custom methods
		void AddCrackTipHeating(void);
		void StartCrackTipHeating(CrackSegment *,Vector &,double);
		
    private:
};

// globals
extern ConductionTask *conduction;

#endif
