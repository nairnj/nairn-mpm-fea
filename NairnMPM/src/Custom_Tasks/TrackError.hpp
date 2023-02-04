/********************************************************************************
    TrackError.hpp
	OSParticulas
    
    Created by Chad Hammerquist on 7/29/2019

    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		CustomTask.hpp
********************************************************************************/

#ifndef _TRACKERRORTASK_

#define _TRACKERRORTASK_

#include "Custom_Tasks/CustomTask.hpp"

#define MAX_COLUMNS 10
class Expression;

class TrackError : public CustomTask
{
    public:
		static int numTETasks;
	
        // constructors and destructors
        TrackError();
        
        // standard methods
		virtual const char *TaskName(void);
		virtual char *InputParam(char *,int &,double &);
		CustomTask * StepCalculation(void);
		virtual CustomTask *Initialize(void);
		virtual CustomTask * PrepareForStep(bool &);
		virtual void SetTextParameter(char *, char *);

	private:
		// settomgs
		int TEnum;
		int numExprs;
		double customArchiveTime, nextCustomArchiveTime;
		int lp_norm_p[MAX_COLUMNS];
		int Cumulative[MAX_COLUMNS];				// 0 false, !=0 is true
		Expression *error_expr[MAX_COLUMNS];

		// run values
		bool doErrorExport;			// flag to export to the file this step
		int count_n;
		char *dummy;
		double total_error[10];
};

#endif

