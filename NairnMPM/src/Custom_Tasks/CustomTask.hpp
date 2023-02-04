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
		static int numberCustomHistoryVariables;
        
        // constructors and destructors
        CustomTask();
        virtual ~CustomTask();
        
        // initialize
		virtual void Finalize(void);
		virtual const char *TaskName(void);
        virtual char *InputParam(char *,int &,double &);
		virtual void SetTextParameter(char *,char *);
        virtual CustomTask *Initialize(void);
		virtual void ClaimHistoryVariables(void);
	
		// run the task
        virtual CustomTask *PrepareForStep(bool &);
        virtual CustomTask *StepCalculation(void);
        virtual CustomTask *Step0Calculation(void);
		virtual CustomTask *FinishForStep(bool &);
        virtual bool HasReport(void);
        virtual CustomTask *Report(void);

		// for particle to node extrapolations
		virtual CustomTask *BeginExtrapolations(void);
		virtual CustomTask *NodalExtrapolation(NodalPoint *,MPMBase *,short,int,double,short);
		virtual CustomTask *EndExtrapolations(void);
    
    protected:
		int ignoreArgument;			// return point for valid parameters with no arguments
};

// globals
extern CustomTask *theTasks;

#endif
