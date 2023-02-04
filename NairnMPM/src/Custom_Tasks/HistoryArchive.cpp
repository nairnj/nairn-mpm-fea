/********************************************************************************
	HistoryArchive.cpp
	nairn-mpm-fea

	Created by John Nairn on 10/26/12.
	Copyright (c) 2012 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Custom_Tasks/HistoryArchive.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "System/ArchiveData.hpp"
#include "System/UnitsController.hpp"

static int historyArg;

#pragma mark Constructors and Destructors

// Constructors
HistoryArchive::HistoryArchive() : CustomTask()
{
	customArchiveTime = -1.;
	nextCustomArchiveTime = -1.;
}

// Return name of this task
const char *HistoryArchive::TaskName(void) { return "Archive particle history data to a tab-delimited file"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
// not thread safe due to push_back()
char *HistoryArchive::InputParam(char *pName,int &input,double &gScaling)
{
    if(strcmp(pName,"archiveTime")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&customArchiveTime,gScaling,1.e-3);
    }
	
    else if(strcmp(pName,"firstArchiveTime")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&nextCustomArchiveTime,gScaling,1.e-3);
    }
	
	else
	{	// assume an integer
		int historyNum;
		sscanf(pName,"%d",&historyNum);
		quantity.push_back(historyNum);
		
		// ignore any command data
		input=INT_NUM;
        return (char *)&historyArg;
	}
	
	// check remaining commands
    return CustomTask::InputParam(pName,input,gScaling);
}

#pragma mark GENERIC TASK METHODS

// called once at start of MPM analysis - initialize and print info
CustomTask *HistoryArchive::Initialize(void)
{
    cout << "Archive particle history results to tab-delimited text files." << endl;
	
	// time interval
	cout << "   Archive time: ";
	if(customArchiveTime>=0.)
	{	cout << customArchiveTime*UnitsController::Scaling(1.e3) << " " << UnitsController::Label(ALTTIME_UNITS);
		if(nextCustomArchiveTime<0.)
		{	nextCustomArchiveTime = customArchiveTime;
			cout << endl;
		}
		else
		{	cout << ", starting at " << nextCustomArchiveTime*UnitsController::Scaling(1.e3) << " " << UnitsController::Label(ALTTIME_UNITS) << endl;
		}
	}
	else
		cout << "same as particle archives" << endl;
	
	// quantities
	unsigned int q;
	cout << "   History Numbers: " ;
	int len=20;
	for(q=0;q<quantity.size();q++)
	{	if(len>70)
		{	cout << "\n      ";
			len=6;
		}
		cout << quantity[q] ;
		if(q < quantity.size()-1) cout << ", " ;
		len += 5;
	}
	cout << endl;
	
    return nextTask;
}

// called when MPM step is getting ready to do custom tasks
// does not use exrapolations so no need to set
CustomTask *HistoryArchive::PrepareForStep(bool &needExtraps)
{
	if(customArchiveTime >= 0.)
	{	if(mtime+timestep >= nextCustomArchiveTime)
		{	doHistoryExport = true;
			nextCustomArchiveTime += customArchiveTime;
		}
		else
			doHistoryExport = false;
	}
	else
		doHistoryExport = archiver->WillArchive();
	if(quantity.size()==0) doHistoryExport = false;
    return nextTask;
}

// Archive VTK file now
CustomTask *HistoryArchive::StepCalculation(void)
{
	if(doHistoryExport)
		archiver->ArchiveHistoryFile(mtime+timestep,quantity);
    return nextTask;
}
