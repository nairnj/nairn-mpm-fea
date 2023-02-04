/********************************************************************************
LoadControl.hpp
nairn-mpm-fea

Created by Chad Hammerquist April 2016
Copyright (c) 2016 John A. Nairn, All rights reserved.

Dependencies
CustomTask.hpp
********************************************************************************/

#ifndef _LOADCONTROLTASK_

#define _LOADCONTROLTASK_

#include "Custom_Tasks/CustomTask.hpp"
class Expression;

class LoadControl : public CustomTask
{
    public:
        // input parameters
        double At,Kp,Ki,Kd;
        double smoothF,smoothA,smoothErr;
        double startTime;
        int material,direction;
        char *matName;
        double velocity,minVelocity;
        Expression *Load_Function;
        char *archive,*archiveBase;

        // constructors and destructors
        LoadControl();

        // standard methods
        virtual const char *TaskName(void);
        virtual char *InputParam(char *, int &, double &);
        virtual CustomTask *Initialize(void);
        virtual CustomTask *StepCalculation(void);
        virtual void SetTextParameter(char *, char *);
    
    private:
        // used in code
        int globalLoadQuantity,globalDistanceQuantity;
        double prevUpdateTime;
        double load,prevLoad;
        int numStartup;
        double sumXY,sumX,sumY,sumX2;
        double error,prevError,intError;
};

#endif
