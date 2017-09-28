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
class MatVelocityField;
class CrackVelocityField;
class NodalValueBC;
class MatPtLoadBC;

enum { CONTACT_EQUILIBRATED=0, CONTACT_CONVECTED };

class TransportTask
{
    public:
		static bool hasContactEnabled;
        TransportTask *nextTask;
       
        // constructors and destructors
        TransportTask();
        virtual ~TransportTask();
        
		// called before first time step, but after preliminary calcs
		virtual TransportTask *Initialize(void) = 0;
	
		// Mass and Momentum extrapolation and post extrapolation for transport property to the grid
		virtual TransportTask *Task1Extrapolation(NodalPoint *,MPMBase *,double,short,int);
		virtual void Task1ContactExtrapolation(NodalPoint *ndpt,short,int,double,double);
		virtual TransportTask *Task1Reduction(NodalPoint *,NodalPoint *);
		virtual void Task1ContactReduction(NodalPoint *,NodalPoint *);
		virtual TransportTask *GetTransportNodalValue(NodalPoint *);
		virtual void GetContactNodalValue(NodalPoint *);
		virtual void ImposeValueBCs(double);
		virtual TransportTask *GetGradients(double);
		virtual void ZeroTransportGradients(MPMBase *) = 0;
		virtual void ZeroTransportContactGradients(MPMBase *);
		virtual void AddTransportGradients(MPMBase *,Vector *,NodalPoint *,short) = 0;
		virtual void AddTransportContactGradients(MPMBase *,Vector *,NodalPoint *,short);
	
		// find forces for transport calculation
		virtual TransportTask *AddForces(NodalPoint *,MPMBase *,double,double,double,double,TransportProperties *,short,int) = 0;
		virtual void AddContactForces(NodalPoint *,MPMBase *,double,double,double,double,TransportProperties *,short,int);
        virtual TransportTask *ForcesReduction(NodalPoint *,NodalPoint *);
		virtual void ForcesContactReduction(NodalPoint *,NodalPoint *);
		virtual TransportTask *SetTransportForceBCs(double);
		virtual void AddFluxCondition(NodalPoint *,double,bool);
		
		// update momentum task and contact flow
		virtual TransportTask *GetTransportRates(NodalPoint *,double);
		virtual void TransportContactRates(NodalPoint *,double);
#ifdef CONTACT_HEAT_FLOW
		// contact calculations
		virtual TransportTask *MatContactFlowCalculations(MatVelocityField *,NodalPoint *,CrackVelocityField *,double,bool);
		virtual TransportTask *CrackContactFlowCalculations(CrackVelocityField *,NodalPoint *,double,int);
#endif
		
		// update particles task
		virtual TransportTask *IncrementTransportRate(NodalPoint *,double,double &,short,int) const;
		virtual TransportTask *MoveTransportValue(MPMBase *,double,double) const;
		
		// update particle strains
		virtual double IncrementValueExtrap(NodalPoint *,double,short,int) const;
		virtual double GetDeltaValue(MPMBase *,double) const;
	
		// if needed for SZS or USAVG, update value on the grid
		virtual TransportTask *UpdateNodalValues(double);
		virtual void UpdateContactNodalValues(NodalPoint *,double);
    
        // accessors
    
        // Return name of this task
        virtual const char *TaskName(void);
        
        // get the next task
        TransportTask *GetNextTransportTask(void) const;
    
        // Get transport mass, mTp in notes, and get particle transport value
        virtual double GetTransportMassAndValue(MPMBase *,double *) = 0;
    
        // pointer to transport field
        virtual TransportField *GetTransportFieldPtr(NodalPoint *) const = 0;
        
        // Get max(kappa)/k for transport task time step calculation
        virtual TransportTask *TransportTimeStepFactor(int,double *) = 0;
    
        // pointer to the first boundary condition
        virtual NodalValueBC *GetFirstBCPtr(void) const = 0;
		virtual MatPtLoadBC *GetFirstFluxBCPtr(void) const = 0;
	
        // get pointers to particle value for this transport property
        virtual double *GetParticleValuePtr(MPMBase *mptr) const = 0;
        virtual double *GetPrevParticleValuePtr(MPMBase *mptr) const = 0;
};

// globals
extern TransportTask *transportTasks;

#endif

