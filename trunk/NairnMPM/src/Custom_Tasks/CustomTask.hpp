/********************************************************************************
    CustomTask.hpp
    NairnMPM
    
    Created by John Nairn on Fri Aug 15 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _CUSTOMTASK_

#define _CUSTOMTASK_

class NodalPoint;
class MPMBase;

class CustomTask
{
    public:
        CustomTask *nextTask;
        
        // constructors and destructors
        CustomTask();
        virtual ~CustomTask();
        
        // standard methods
		virtual const char *TaskName(void);
        virtual char *InputParam(char *,int &);
        virtual CustomTask *Initialize(void);
        virtual CustomTask *PrepareForStep(bool &);
        virtual CustomTask *FinishForStep(void);
        virtual CustomTask *BeginExtrapolations(void);
        virtual CustomTask *EndExtrapolations(void);
        virtual CustomTask *NodalExtrapolation(NodalPoint *,MPMBase *,short,int,double);
        virtual CustomTask *ParticleCalculation(NodalPoint *,MPMBase *,short,int,double,double,double,double);
        virtual CustomTask *ParticleExtrapolation(MPMBase *);
        virtual CustomTask *StepCalculation(void);
        
    private:
};

// globals
extern CustomTask *theTasks;

#endif
