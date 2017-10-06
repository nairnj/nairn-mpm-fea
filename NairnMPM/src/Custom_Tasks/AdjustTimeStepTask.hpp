/********************************************************************************
    AdjustTimeStepTask.hpp
    nairn-mpm-fea

    Created by John Nairn on 9/24/12.
    Copyright (c) 2012 John A. Nairn, All rights reserved.

    Dependencies
        CustomTask.hpp
********************************************************************************/

#ifndef _ADJUSTTIMESTEPTASK_

#define _ADJUSTTIMESTEPTASK_

#include "Custom_Tasks/CustomTask.hpp"

class AdjustTimeStepTask : public CustomTask
{
    public:
    
        // constructors and destructors
        AdjustTimeStepTask();
    
        // standard methods
        virtual const char *TaskName(void);
        virtual char *InputParam(char *,int &,double &);
        virtual CustomTask *Initialize(void);
	
        virtual CustomTask *PrepareForStep(bool &);
        virtual CustomTask *StepCalculation(void);
    
    private:
        double customAdjustTime,nextCustomAdjustTime;
        bool doAdjust;			// flag to adjust time step this step
        int verbose;            // 0 or 1 to report changes
        double lastReportedTimeStep;
		double velocityCFL;
    	double reportRatio;
    
};

#endif

