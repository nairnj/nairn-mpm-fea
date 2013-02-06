/********************************************************************************
    AdjustTimeStepTask.cpp
    NairnMPM

    Created by John Nairn on 9/24/12.
    Copyright (c) 2012 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Custom_Tasks/AdjustTimeStepTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "System/ArchiveData.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "Elements/ElementBase.hpp"

#pragma mark Constructors and Destructors

// Constructors
AdjustTimeStepTask::AdjustTimeStepTask()
{
	customAdjustTime = -1.;
	nextCustomAdjustTime = -1.;
    verbose = 0;
    lastReportedTimeStep = -1;
}

// Return name of this task
const char *AdjustTimeStepTask::TaskName(void) { return "Periodically adjust MPM time step"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
char *AdjustTimeStepTask::InputParam(char *pName,int &input)
{
    if(strcmp(pName,"adjustTime")==0)
    {	input=DOUBLE_NUM;
        return (char *)&customAdjustTime;				// assumes in ms
    }
		
    if(strcmp(pName,"verbose")==0)
    {	input=INT_NUM;
        return (char *)&verbose;
    }
	// check remaining commands
    return CustomTask::InputParam(pName,input);
}

#pragma mark GENERIC TASK METHODS

// called once at start of MPM analysis - initialize and print info
CustomTask *AdjustTimeStepTask::Initialize(void)
{
    cout << "Periodically adjust MPM time step." << endl;
	
	// time interval
	cout << "   Adjust interval: ";
	if(customAdjustTime>0)
	{	cout << customAdjustTime << " ms" << endl;
		customAdjustTime/=1000.;				// convert to sec
		nextCustomAdjustTime=customAdjustTime;
	}
	else
		cout << "same as particle archives" << endl;
    if(verbose!=0)
        cout << "   Verbose: yes" << endl;
    else
        cout << "   Verbose: no" << endl;
	
    return nextTask;
}

// called when MPM step is getting ready to do custom tasks
CustomTask *AdjustTimeStepTask::PrepareForStep(bool &needExtraps)
{
	if(customAdjustTime>0.)
	{	if(mtime+timestep>=nextCustomAdjustTime)
        {	doAdjust=TRUE;
            nextCustomAdjustTime+=customAdjustTime;
        }
        else
            doAdjust=FALSE;
	}
	else
		doAdjust=archiver->WillArchive();
    
    needExtraps = FALSE;
    return nextTask;
}

// Adjust time step now
CustomTask *AdjustTimeStepTask::StepCalculation(void)
{
    // exit when not needed
    if(!doAdjust) return nextTask;
    
    int p;
    short matid;
    double area,volume,crot,tst,dcell;
    
    // get grid dimensions
    double minSize=mpmgrid.GetMinCellDimension()/10.;	// in cm
    
    if(lastReportedTimeStep<0) lastReportedTimeStep = timestep;
    
    // reset globals
    timestep=1.e15;
    propTime=1.e15;
    
    // loop over material points
    for(p=0;p<nmpms;p++)
	{	// verify material is defined and sets if field number (in in multimaterial mode)
		matid=mpm[p]->MatID();
		
		// skip rigid materials
		if(theMaterials[matid]->Rigid()) continue;
        
		// element and mp properties
		if(fmobj->IsThreeD())
		{	volume=theElements[mpm[p]->ElemID()]->GetVolume()/1000.;	// in cm^3
			dcell = (minSize>0.) ? minSize : pow(volume,1./3.) ;
		}
		else
		{	area=theElements[mpm[p]->ElemID()]->GetArea()/100.;	// in cm^2
			volume=mpm[p]->thickness()*area/10.;				// in cm^2
			dcell = (minSize>0.) ? minSize : sqrt(area) ;
		}
        
        // check time step using convergence condition
        crot=theMaterials[matid]->CurrentWaveSpeed(fmobj->IsThreeD(),mpm[p])/10.;           // in cm/sec
		tst=fmobj->GetCFLCondition()*dcell/crot;                                            // in sec
        if(tst<timestep) timestep=tst;
        
        // propagation time (in sec)
        tst=fmobj->GetPropagationCFLCondition()*dcell/crot;                                        // in sec
        if(tst<propTime) propTime=tst;
	}
	
    // verify time step and make smaller if needed
	strainTimestep = (fmobj->mpmApproach==USAVG_METHOD) ? timestep/2. : timestep ;
    
	// propagation time step (no less than timestep)
    if(propTime<timestep) propTime=timestep;
    
    // report if changed by 5% since last reported change
    if(verbose!=0)
    {   double ratio = timestep/lastReportedTimeStep;
        if(ratio < 0.95)
            cout << "# time step reduced to " << timestep*1000. << " ms" << endl;
        else if(ratio > 1.05)
            cout << "# time step increased to " << timestep*1000. << " ms" << endl;
        lastReportedTimeStep = timestep;
    }

    return nextTask;
}

