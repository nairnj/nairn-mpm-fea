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

// To extrpolate c/crel for materials where crel changes with position (in gTValueRel)
// If comment out, extrapolates just crel (in gTRelValue) and gets gradient from
//		extrapolated c/extrapolated crel
#define USE_GTVALUEREL

// to allow c/csat>1, but not Less than zero (only affects diffusion code)
//#define NO_UPPER_LIMIT

enum { CONTACT_EQUILIBRATED=0, CONTACT_CONVECTED };

// One for conduction, one for diffusion (or poroelastity), and
// one for each other task available in code (currently phase field diffusion)
#define MAX_TRANSPORT_TASKS 3
#define MAX_DIFFUSION_TASKS 2

class TransportTask
{
    public:
		static bool hasContactEnabled;
        TransportTask *nextTask;
		double transportTimeStep;
		static int XPICOrder;
		static bool hasXPICOption;
		bool usingXPIC;
		double usingFraction;
		bool doCrelExtrapolation;
       
        // constructors and destructors
        TransportTask();
        virtual ~TransportTask();
		void CheckTimeStep(double);
		double GetTimeStep(void) const;
	
		// called before first time step, but after preliminary calcs
		virtual TransportTask *Initialize(void) = 0;
	
		// Mass and Momentum extrapolation and post extrapolation for transport property to the grid
		virtual TransportTask *Task1Extrapolation(NodalPoint *,MPMBase *,double,short,int) = 0;
		virtual double GetVpCTp(MPMBase *) = 0;
		virtual void Task1ContactExtrapolation(NodalPoint *ndpt,short,int,double,double);
		virtual TransportTask *Task1Reduction(NodalPoint *,NodalPoint *);
		virtual void Task1ContactReduction(NodalPoint *,NodalPoint *);
		virtual TransportTask *GetTransportNodalValue(NodalPoint *);
		virtual void GetContactNodalValue(NodalPoint *);
		virtual void ImposeValueBCs(double,bool);
#ifdef TRANSPORT_FMPM
		virtual TransportTask *RestoreValueBCs(void);
		virtual TransportTask *ImposeValueGridBCs(double,double,int);
		virtual TransportTask *SetTransportFluxBCs(void);
#else
		virtual TransportTask *SetTransportForceAndFluxBCs(double);
#endif
		virtual TransportTask *GetGradients(double);
		virtual void ZeroTransportGradients(MPMBase *) = 0;
		virtual void AddTransportGradients(MPMBase *,Vector *,NodalPoint *,short) = 0;
		virtual void ZeroTransportContactGradients(MPMBase *);
		virtual void AddTransportContactGradients(MPMBase *,Vector *,NodalPoint *,short);
	
		// find forces for transport calculation
		virtual TransportTask *AddForces(NodalPoint *,MPMBase *,double,double,double,double,TransportProperties *,short,int) = 0;
		virtual void AddContactForces(NodalPoint *,MPMBase *,double,double,double,double,TransportProperties *,short,int);
        virtual TransportTask *ForcesReduction(NodalPoint *,NodalPoint *);
		virtual void ForcesContactReduction(NodalPoint *,NodalPoint *);
		virtual void AddFluxCondition(NodalPoint *,double,bool);
		
		// update momentum task and contact flow
		virtual TransportTask *UpdateTransport(NodalPoint *,double);
		virtual void TransportContactRates(NodalPoint *,double);

		// update particles task
		virtual TransportTask *InitializeForXPIC(NodalPoint *,double,int) const;
		virtual double IncrementTransportRate(NodalPoint *,double,short,int) const;
		virtual TransportTask *MoveTransportValue(MPMBase *,double,double,double) const;
		virtual void AdjustRateAndValue(MPMBase *,double &,double &,double &,double) const;

		// update particle strains
		virtual double IncrementValueExtrap(NodalPoint *,double,short,int) const;
		virtual double IncrementLumpedValueExtrap(NodalPoint *,double,short,int) const;
		virtual void GetDeltaValue(MPMBase *,double,double *) const;
	
        // accessors
    
        // Return name of this task
        virtual const char *TaskName(void);
        
        // get the next task
        TransportTask *GetNextTransportTask(void) const;
    
        // pointer to transport field
        virtual TransportField *GetTransportFieldPtr(NodalPoint *) const = 0;
        
        // Get max(kappa)/k for transport task time step calculation
        virtual TransportTask *TransportTimeStepFactor(int,double *) = 0;
    
        // pointer to the first boundary condition
        virtual NodalValueBC *GetFirstBCPtr(void) const = 0;
		virtual MatPtLoadBC *GetFirstFluxBCPtr(void) const = 0;
	
        // get pointers to particle value for this transport property
		virtual double GetParticleValue(MPMBase *mptr) const = 0;
        virtual double *GetParticleValuePtr(MPMBase *mptr) const = 0;
		virtual double GetPrevParticleValue(MPMBase *mptr) const = 0;
        virtual double *GetPrevParticleValuePtr(MPMBase *mptr) const = 0;
	
		// if using XPIC
		void SetUsingTransportXPIC(bool,double);
		virtual bool IsUsingTransportXPIC(void) const;
		virtual bool IsUsingTransportXPIC(double &) const;
		virtual bool ShouldBlendFromGrid(double &) const;

        // multiple tasks of same type (diffusion only, but here to allow access)
        virtual int GetNumber(void) const;
        virtual void SetNumber(int);
    
		// static methods
		static void GetTransportValues(NodalPoint *);
#ifdef TRANSPORT_FMPM
		static void TransportGridBCs(double,double,int);
#endif
		static void TransportBCsAndGradients(double);
		static void UpdateTransportOnGrid(NodalPoint *);
		static void TransportForceBCs(double);

};

// globals
extern TransportTask *transportTasks;
extern int numTransport;


#endif

