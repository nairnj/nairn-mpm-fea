/********************************************************************************
    CustomTask.hpp
    nairn-mpm-fea
    
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
        
        // initialize
		virtual const char *TaskName(void);
        virtual char *InputParam(char *,int &);
        virtual CustomTask *Initialize(void);
	
		// run the task
        virtual CustomTask *PrepareForStep(bool &);
        virtual CustomTask *StepCalculation(void);
		virtual CustomTask *FinishForStep(void);
	
		// for particle to node extrapolations
		virtual CustomTask *BeginExtrapolations(void);
		virtual CustomTask *NodalExtrapolation(NodalPoint *,MPMBase *,short,int,double,short);
		virtual CustomTask *EndExtrapolations(void);
    
    private:
};

// globals
extern CustomTask *theTasks;

#endif
