/********************************************************************************
	SetCustomTasks.cpp - extra code for MPMReadHandler.cpp for creating custom tasks
	nairn-mpm-fea

	Created by John Nairn on July 19, 2018.
	Copyright (c) 2018 RSAC Software. All rights reserved.
 ********************************************************************************/

#include "stdafx.h"
#include "Read_MPM/MPMReadHandler.hpp"
#include "Custom_Tasks/ConductionTask.hpp"

// Custom tasks
#include "Custom_Tasks/ReverseLoad.hpp"
#include "Custom_Tasks/VTKArchive.hpp"
#include "Custom_Tasks/HistoryArchive.hpp"
#include "Custom_Tasks/AdjustTimeStepTask.hpp"
#include "Custom_Tasks/CarnotCycle.hpp"
#include "Custom_Tasks/CustomThermalRamp.hpp"
#include "Custom_Tasks/PeriodicXPIC.hpp"
#include "Custom_Tasks/DeleteDamaged.hpp"
#include "Custom_Tasks/FluidSource.hpp"
#include "Custom_Tasks/TrackError.hpp"
#include "Custom_Tasks/LoadControl.hpp"

// Create custom task
void MPMReadHandler::ScheduleCustomTask(const Attributes& attrs)
{
	char *aName,*value;
	unsigned int i,numAttr;
	
	CustomTask *nextTask=NULL;
	numAttr = (int)attrs.getLength();
	
	for(i=0;i<numAttr;i++)
	{   aName=XMLString::transcode(attrs.getLocalName(i));
		if(strcmp(aName,"name")==0)
		{   value=XMLString::transcode(attrs.getValue(i));
			
			if(strcmp(value,"ReverseLoad")==0)
			{   nextTask=(CustomTask *)(new ReverseLoad());
				if(nextTask==NULL) throw SAXException("Out of memory creating a custom task.");
			}
			
			else if(strcmp(value,"VTKArchive")==0)
			{   nextTask=(CustomTask *)(new VTKArchive());
				if(nextTask==NULL) throw SAXException("Out of memory creating a custom task.");
			}
			
			else if(strcmp(value,"HistoryArchive")==0)
			{   nextTask=(CustomTask *)(new HistoryArchive());
				if(nextTask==NULL) throw SAXException("Out of memory creating a custom task.");
			}
			
			else if(strcmp(value,"AdjustTimeStep")==0)
			{   nextTask=(CustomTask *)(new AdjustTimeStepTask());
				if(nextTask==NULL) throw SAXException("Out of memory creating a custom task.");
			}
			
			else if(strcmp(value,"CarnotCycle")==0)
			{   nextTask=(CustomTask *)(new CarnotCycle());
				if(nextTask==NULL) throw SAXException("Out of memory creating a custom task.");
			}
			
			else if(strcmp(value,"PropertyRamp")==0 || strcmp(value,"ThermalRamp")==0)
			{   nextTask=(CustomTask *)(new CustomThermalRamp());
				if(nextTask==NULL) throw SAXException("Out of memory creating a custom task.");
				ConductionTask::activeRamp = true;
			}

			else if(strcmp(value,"PeriodicXPIC")==0)
			{   nextTask=(CustomTask *)(new PeriodicXPIC());
				if(nextTask==NULL) throw SAXException("Out of memory creating a custom task.");
			}
			
            else if (strcmp(value, "DeleteDamaged") == 0)
            {    nextTask = (CustomTask *)(new DeleteDamaged());
                if (nextTask == NULL) throw SAXException("Out of memory creating a custom task.");
            }
            
            else if (strcmp(value, "FluidSource") == 0)
            {    nextTask = (CustomTask *)(new FluidSource());
                if (nextTask == NULL) throw SAXException("Out of memory creating a custom task.");
            }
            
            else if (strcmp(value, "TrackError") == 0)
            {   nextTask = (CustomTask *)(new TrackError());
                if (nextTask == NULL) throw SAXException("Out of memory creating a custom task.");
            }
            
            else if(strcmp(value,"LoadControl")==0)
            {   nextTask=(CustomTask *)(new LoadControl());
                if(nextTask==NULL) throw SAXException("Out of memory creating a custom task.");
            }
            
			else
				throw SAXException("Unknown custom task requested for scheduling.");
			
			delete [] value;
		}
		delete [] aName;
		
		// exit when find a match
		if(nextTask!=NULL) break;
	}
	
	// link into tasks list
	if(nextTask != NULL)
	{	// set block to all taks parameters
		block = TASKPARAMETERS;
		
		// start list or add after previous custom task
		if (currentTask == NULL)
		{	theTasks = nextTask;
		}
		else
		{	//cout << ((CustomTask *)currentTask)->TaskName()<<endl;
			((CustomTask *)currentTask)->nextTask = nextTask;
			//cout << nextTask->TaskName() << endl;
		}
		
		// set current task
		currentTask=(char *)nextTask;
	}
}

// set custom task parameter
void MPMReadHandler::SetCustomTasksParameter(const Attributes& attrs)
{
	char *aName,*value;
	unsigned int i,numAttr;
	
	inputPtr=NULL;
	numAttr = (int)attrs.getLength();
	for(i=0;i<numAttr;i++)
	{   aName=XMLString::transcode(attrs.getLocalName(i));
		if(strcmp(aName,"name")==0)
		{	// check for parameter
			value=XMLString::transcode(attrs.getValue(i));
			inputPtr=((CustomTask *)currentTask)->InputParam(value,input,gScaling);
			
			// error message
			if(inputPtr==NULL)
			{	char errMsg[200];
                size_t errSize=200;
				snprintf(errMsg,errSize,"Unrecognized custom task parameter: '%s'",value);
				delete [] value;
				throw SAXException(errMsg);
			}
				
			delete [] value;
		}
		delete [] aName;
		if(inputPtr!=NULL) break;
	}
	
	// if still NULL, had no name attribute
	if(inputPtr==NULL)
		throw SAXException("Custom task parameter has no name.");
}

