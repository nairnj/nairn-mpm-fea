/********************************************************************************
	CustomThermalRamp.hpp
	nairn-mpm-fea

	Created by John Nairn on 9/14/16.
	Copyright (c) 2016 John A. Nairn, All rights reserved.

	Dependencies
		CustomTask.hpp
********************************************************************************/

#ifndef _THERMALRAMPTASK_

#define _THERMALRAMPTASK_

#include "Custom_Tasks/CustomTask.hpp"
class ROperation;
class BMPLevel;

class CustomThermalRamp : public CustomTask
{
	public:
	
		// constructors and destructors
		CustomThermalRamp();
		virtual ~CustomThermalRamp();
	
		// standard methods
		virtual const char *TaskName(void);
		virtual char *InputParam(char *,int &,double &);
		virtual void SetTextParameter(char *,char *);
		virtual CustomTask *Initialize(void);
	
		virtual CustomTask *PrepareForStep(bool &);
		virtual CustomTask *StepCalculation(void);
		virtual CustomTask *FinishForStep(bool &);
	
	private:
		double isoDeltaT;		// final temperature change
		double isoRampTime;		// time to reach temperature change
		double rampStart;		// when to start the ramp
		double endTime;			// all done by this time
		double currentDeltaT;	// track Delta T to allow variable time steps
		bool doRamp;			// flag to adjust temperature this time step
		int sigmoidal;			// Use sigmoidal shape
		ROperation *scaleFxn;
		static double varTime,varXValue,varYValue,varZValue;
	
		// for imaged ramp
		char *bmpFile;
		Vector orig;
		double width,height;
		bool flipped;
		double deltaTmin,deltaTmax,deltaTscale;
		BMPLevel *firstLevel;
		BMPLevel *currentLevel;
		BMPLevel *firstMatPt;
		BMPLevel *currentMatPt;
};

#endif

