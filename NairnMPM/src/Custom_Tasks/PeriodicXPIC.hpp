/********************************************************************************
	Periodic XPIC custom Custon Task
	nairn-mpm-fea
 
	Created by John Nairn on 1/26/18.
	Copyright (c) 2018 John A. Nairn, All rights reserved.
 
	Dependencies
		CustomTask.hpp
********************************************************************************/

#ifndef _PERIODICXPICTASK_

#define _PERIODICXPICTASK_

#include "Custom_Tasks/CustomTask.hpp"
#include "Custom_Tasks/TransportTask.hpp"

class PeriodicXPIC : public CustomTask
{
	public:
	
		// constructors and destructors
		PeriodicXPIC();
	
		// standard methods
		virtual void Finalize(void);
		virtual const char *TaskName(void);
		virtual char *InputParam(char *,int &,double &);
		virtual void SetTextParameter(char *,char *);
		virtual CustomTask *Initialize(void);
	
		virtual CustomTask *PrepareForStep(bool &);
		virtual CustomTask *StepCalculation(void);
	
		virtual const char *GetType(void);
	
	private:
		// parameters
		double periodicTime;
		double periodicCFL;
		int periodicSteps,periodicXPICorder;
		int verbose;
		double nextPeriodicTime;
		int nextPeriodicStep;
		bool doXPIC;
		bool usingFMPM;
		int gridBCOption,periodicFMPMorder;
		double periodicTimeConduction,periodicTimeDiffusion[MAX_DIFFUSION_TASKS];
		double periodicCFLConduction,periodicCFLDiffusion[MAX_DIFFUSION_TASKS];
		int periodicStepsConduction,periodicStepsDiffusion[MAX_DIFFUSION_TASKS];;
		double conductionFractionFMPM,diffusionFractionFMPM[MAX_DIFFUSION_TASKS];
		double nextPeriodicTimeConduction,nextPeriodicTimeDiffusion[MAX_DIFFUSION_TASKS];
		int nextPeriodicStepConduction,nextPeriodicStepDiffusion[MAX_DIFFUSION_TASKS];
		bool doXPICConduction,doXPICDiffusion[MAX_DIFFUSION_TASKS];
};

#endif

