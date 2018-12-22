/********************************************************************************
    ConductionTask.hpp
    nairn-mpm-fea
    
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
class MatVelocityField;
class NodalValueBC;

class ConductionTask : public TransportTask
{
    public:
		static bool active,activeRamp,crackTipHeating,adiabatic;
		static bool crackContactHeating,matContactHeating;
	
		// initialize
		virtual TransportTask *Initialize(void);
	
		// mass and momentum
		virtual TransportTask *Task1Extrapolation(NodalPoint *,MPMBase *,double,short,int);
		virtual double GetVpCTp(MPMBase *);
		virtual void ZeroTransportGradients(MPMBase *);
		virtual void AddTransportGradients(MPMBase *,Vector *,NodalPoint *,short);
	
		// grid force extrapolation
		virtual TransportTask *AddForces(NodalPoint *,MPMBase *,double,double,double,double,TransportProperties *,short,int);
		virtual TransportTask *SetTransportForceBCs(double);
	
		// custom methods
		void AddCrackTipHeating(void);
		void StartCrackTipHeating(CrackSegment *,Vector &,double);
    
        // accessors
        virtual const char *TaskName(void);
        virtual TransportTask *TransportTimeStepFactor(int,double *);
        virtual TransportField *GetTransportFieldPtr(NodalPoint *) const;
        virtual NodalValueBC *GetFirstBCPtr(void) const;
		virtual MatPtLoadBC *GetFirstFluxBCPtr(void) const;
        virtual double *GetParticleValuePtr(MPMBase *mptr) const;
        virtual double *GetPrevParticleValuePtr(MPMBase *mptr) const;
    
        // class methods
        static void ThermodynamicsOutput(void);
        static bool IsSystemIsolated(void);
		
    protected:
		// contact heat flow
		int crackGradT,materialGradT;
};

// globals
extern ConductionTask *conduction;

#endif
