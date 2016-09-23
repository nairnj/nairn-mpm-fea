/********************************************************************************
	CustomThermalRamp.cpp
	nairn-mpm-fea

	Created by John Nairn on 9/14/16.
	Copyright (c) 2016 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Custom_Tasks/CustomThermalRamp.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "System/UnitsController.hpp"
#include "Read_XML/mathexpr.hpp"

extern double timestep;

double CustomThermalRamp::varTime=0.;
double CustomThermalRamp::varXValue=0.;
double CustomThermalRamp::varYValue=0.;
double CustomThermalRamp::varZValue=0.;
PRVar rampVarArray[4] = { NULL, NULL, NULL, NULL };

#pragma mark Constructors and Destructors

// Constructors
CustomThermalRamp::CustomThermalRamp()
{
	isoRampTime = -1.;
	rampStart = 0.;
	isoDeltaT = 0.;
	sigmoidal = 0;
	scaleFxn = NULL;
}

// Destructor (and it is virtual)
CustomThermalRamp::~CustomThermalRamp()
{	if(scaleFxn!=NULL)
	{	delete rampVarArray[0];
		delete rampVarArray[1];
		delete rampVarArray[2];
		delete rampVarArray[3];
		delete scaleFxn;
	}
}

// Return name of this task
const char *CustomThermalRamp::TaskName(void) { return "Ramp particle temperatures"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
char *CustomThermalRamp::InputParam(char *pName,int &input,double &gScaling)
{
	if(strcmp(pName,"time")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&isoRampTime,gScaling,1.e-3);
	}
	
	else if(strcmp(pName,"DeltaT")==0)
	{	input=DOUBLE_NUM;
		return (char *)&isoDeltaT;
	}
	
	else if(strcmp(pName,"start")==0)
	{	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&rampStart,gScaling,1.e-3);
	}
	
	else if(strcmp(pName,"sigmoidal")==0)
	{	input=INT_NUM;
		return (char *)&sigmoidal;
	}
	
	else
	{	char* check4function = strstr(pName,"scale");
		if(check4function==&pName[0])
		{	char *expr = new char[strlen(pName)];
			if(strlen(pName)<7)
				strcpy(expr,"");
			else
				strcpy(expr,&pName[6]);
			SetTextParameter(expr); // actually sets function
			delete [] expr;
			input=INT_NUM;
			return (char *)&ignoreArgument;
		}
	}
	
	// check remaining commands
	return CustomTask::InputParam(pName,input,gScaling);
}

// Actually sets load function
// throws std::bad_alloc, SAXException()
void CustomThermalRamp::SetTextParameter(char *fxn)
{
	if(scaleFxn!=NULL)
		ThrowSAXException("Duplicate ramp scale function was entered");
	if(fxn==NULL)
		ThrowSAXException("Ramp scale function is missing");
	if(strlen(fxn)==0)
		ThrowSAXException("Ramp scale function is missing");
	
	// create variable
	if(rampVarArray[0]==NULL)
	{	rampVarArray[0]=new RVar("t",&varTime);
		rampVarArray[1]=new RVar("x",&varXValue);
		rampVarArray[2]=new RVar("y",&varYValue);
		rampVarArray[3]=new RVar("z",&varZValue);
	}
	
	// create function
	scaleFxn = new ROperation(fxn,4,rampVarArray);
	if(scaleFxn->HasError())
		ThrowSAXException("Ramp scale function is not valid");
}

#pragma mark GENERIC TASK METHODS

// called once at start of MPM analysis - initialize and print info
CustomTask *CustomThermalRamp::Initialize(void)
{
	// Ramp should be at least one time step
	int nsteps = (int)(isoRampTime/timestep);
	if(nsteps<1)
	{	nsteps = 1;
		isoRampTime = timestep;
	}
	
	// start after zero
	if(rampStart<0.) rampStart = 0.;
	
	// official end time
	endTime = rampStart+isoRampTime+timestep;
	
	// current Delta T
	currentDeltaT = 0.;
	
	cout << "Ramp particle temperatures." << endl;
	
	char hline[200];
	sprintf(hline,"   Final isothermal temperature difference: %g C",isoDeltaT);
	cout << hline << endl;
	sprintf(hline,"   Ramped between %g and %g %s",rampStart*UnitsController::Scaling(1.e3),
				(rampStart+isoRampTime)*UnitsController::Scaling(1.e3),UnitsController::Label(ALTTIME_UNITS));
	cout << hline << endl;
	cout << "      (which covers " << nsteps << " time steps)" << endl;
	if(sigmoidal)
		cout << "      (use sigmoidal ramp)" << endl;
	if(scaleFxn!=NULL)
	{	char *expr = scaleFxn->Expr('#');
		cout << "   Scaling function =  " << expr << endl;
		delete [] expr;
	}
	
	return nextTask;
}

// called when MPM step is getting ready to do custom tasks
// never uses extrapolations so no need to set needExtraps
CustomTask *CustomThermalRamp::PrepareForStep(bool &needExtraps)
{
	double effTime = mtime+timestep;
	if(effTime<rampStart || effTime>endTime)
		doRamp = false;
	else
		doRamp = true;
	
	return nextTask;
}

// Adjust time step now
CustomTask *CustomThermalRamp::StepCalculation(void)
{
	// exit when not needed
	if(!doRamp) return nextTask;
	
	// get new temperature
	double rampFraction = (mtime+timestep-rampStart)/isoRampTime;
	if(rampFraction>1.)
	{	rampFraction = 1.;
	}
	else if(sigmoidal)
	{	rampFraction = 1./(1+exp(-12.*(rampFraction-0.5)));
	}
	double newDeltaT = rampFraction*isoDeltaT;
	
	// get temperature change
	double deltaT = newDeltaT - currentDeltaT;
	currentDeltaT = newDeltaT;
	
	// loop over nonrigid material points and update temperature
	for(int p=0;p<nmpmsNR;p++)
	{	if(scaleFxn==NULL)
			mpm[p]->pTemperature += deltaT;
		else
		{	varTime = mtime*UnitsController::Scaling(1000.);
			varXValue = mpm[p]->pos.x;
			varYValue = mpm[p]->pos.y;
			varZValue = mpm[p]->pos.y;
			mpm[p]->pTemperature += deltaT*scaleFxn->Val();
		}
	}
	
	return nextTask;
}

// Called when custom tasks are all done on a step
CustomTask *CustomThermalRamp::FinishForStep(bool &removeMe)
{	double effTime = mtime+timestep;
	if(effTime>endTime) removeMe = true;
	return nextTask;
}


