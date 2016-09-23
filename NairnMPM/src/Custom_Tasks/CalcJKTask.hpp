/********************************************************************************
    CalcJKTask.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Fri Aug 15 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.

	Dependencies
		CustomTask.hpp
********************************************************************************/

#ifndef _CALCJKTASK_

#define _CALCJKTASK_

#include "Custom_Tasks/CustomTask.hpp"

class CalcJKTask : public CustomTask
{
    public:
        
        // constructors and destructors
        CalcJKTask();
        
        // standard methods
		virtual const char *TaskName(void);
        virtual CustomTask *Initialize(void);
	
        virtual CustomTask *PrepareForStep(bool &);
        virtual CustomTask *StepCalculation(void);
		virtual CustomTask *FinishForStep(bool &);
    
        // other methods
        void ScheduleJK(int);

    private:
        int getJKThisStep;
};

extern CalcJKTask *theJKTask;

#endif
