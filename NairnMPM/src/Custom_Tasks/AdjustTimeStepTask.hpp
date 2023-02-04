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
    
        virtual void SetLastReportedTimeStep(double);
        virtual bool HasReport(void);
        virtual CustomTask *Report(void);

        // class methods
        static void ChangeTimestep(double,double,bool);
    
    private:
        double customAdjustTime,nextCustomAdjustTime;
        bool doAdjust;			// flag to adjust time step this step
        int verbose;            // 0 or 1 to report changes
        int checkTransportTimeStep;   // 0 or 1 to check transport time step too
        double lastReportedTimeStep;
		double velocityCFL;
    	double reportRatio;
        double maxIncrease;
    
};

extern AdjustTimeStepTask *adjustTimeStepTask;


#endif

