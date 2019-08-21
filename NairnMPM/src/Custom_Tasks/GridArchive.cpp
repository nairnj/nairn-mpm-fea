/********************************************************************************
	GridArchive.cpp
	NairnMPM

	Created by John Nairn on 11/4/14.
	Parent class for custom tasks to archive grid results
********************************************************************************/

#include "stdafx.h"
#include "Custom_Tasks/GridArchive.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "System/ArchiveData.hpp"
#include "System/UnitsController.hpp"

#pragma mark Constructors and Destructors

// Constructors
GridArchive::GridArchive()
{
	customArchiveTime = -1.;          // input in ms, stored in sec
	nextCustomArchiveTime = -1.;      // input in ms, stored in sec
	thisMaterial = -1;				  // all by default
}

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
// not thread safe due to push_back()
char *GridArchive::InputParam(char *pName,int &input,double &gScaling)
{
    if(strcmp(pName,"archiveTime")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&customArchiveTime,gScaling,1.e-3);
    }
	
    else if(strcmp(pName,"firstArchiveTime")==0)
    {	input=DOUBLE_NUM;
		return UnitsController::ScaledPtr((char *)&nextCustomArchiveTime,gScaling,1.e-3);
    }
	
	// trap VTK archive parameters and return if found
	else if(strcmp(pName, "selectMaterial") == 0)
	{	input = INT_NUM;
		return (char *)&thisMaterial;
	}
	
	// check remaining commands in super class
    return CustomTask::InputParam(pName,input,gScaling);
}

#pragma mark GENERIC TASK METHODS

// called once at start of MPM analysis - initialize and print info
// This parent class reports archiving intervals and retures the next task
CustomTask *GridArchive::Initialize(void)
{
	// time interval
	cout << "   Archive time: ";
	if(customArchiveTime>=0.)
	{	cout << customArchiveTime*UnitsController::Scaling(1.e3) << " " << UnitsController::Label(ALTTIME_UNITS);
		if(nextCustomArchiveTime<0.)
		{	nextCustomArchiveTime = 0.0;
			cout << endl;
		}
		else
		{	cout << ", starting at " << nextCustomArchiveTime*UnitsController::Scaling(1.e3) << " " << UnitsController::Label(ALTTIME_UNITS) << endl;
		}
	}
	else
		cout << "same as particle archives" << endl;
	
	
    return nextTask;
}

// called when MPM step is getting ready to do custom tasks
CustomTask *GridArchive::PrepareForStep(bool &needExtraps)
{
    // see if need to export on this time step
	if(customArchiveTime>=0.)
	{	if(mtime+timestep>=nextCustomArchiveTime)
		{	doExport = true;
			nextCustomArchiveTime += customArchiveTime;
		}
		else
			doExport=FALSE;
	}
    else if(mtime<0.5*timestep)
        doExport = TRUE;
	else
		doExport=archiver->WillArchive();
	
    // if no buffer, then no need to do nodal extrapolations
	getExtrapolations = CheckExportForExtrapolations();
	if(getExtrapolations) needExtraps = true;
	
    // retrn next task
    return nextTask;
}

// It is called once (during PrepareForStep()) and doExport will be true or false is scheduled to archive.
//	  A subclass can cancel or force archiving by changind doExport.
// Return true or false if export will need extrapolations (but false if doExport is false)
bool GridArchive::CheckExportForExtrapolations(void) { return doExport; }

// Archive VTK file now
CustomTask *GridArchive::StepCalculation(void)
{	// exit if not export this time step
	if(!doExport) return nextTask;
	
	// tell subclass to export
	ExportExtrapolationsToFiles();
	
	// return nextTasks when odne
    return nextTask;
}

// The extrapolations are done and can no be written to a file
void GridArchive::ExportExtrapolationsToFiles(void) {}

#pragma mark TASK EXTRAPOLATION METHODS

// initialize buffuers (if needed) for grid extrapolations
CustomTask *GridArchive::BeginExtrapolations(void)
{	if(!getExtrapolations) return nextTask;
	AllocateExtrapolationBuffers();
	return nextTask;
}

// subclasses will need to override and allocate buffers needed
// for extrapolation to the grid
void GridArchive::AllocateExtrapolationBuffers(void) {}

// Extrapolate particle data for particle *mpnt to nodal point *ndmi
// vfld and matfld are the crack velocity and material velocity field for this particle-node paiur
// wt is extraolation weight and equal to mp*Sip (or particle mass times the shape function)
// isRigid will be true or false if material for this particle is a rigid contact particle
//CustomTask *VTKArchive::NodalExtrapolation(NodalPoint *ndmi,MPMBase *mpnt,short vfld,int matfld,double wt,short isRigid)

// Finish up extrapolation calculations
CustomTask *GridArchive::EndExtrapolations(void)
{	if(!getExtrapolations) return nextTask;
	FinishExtrapolationCalculations();
	return nextTask;
}

// When extrapolations are done, do any remaining calculations, such as to divide by nodal mass
void GridArchive::FinishExtrapolationCalculations(void) {}


