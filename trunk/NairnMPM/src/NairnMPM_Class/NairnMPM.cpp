/*********************************************************************
    NairnMPM.cpp
    Nairn Research Group MPM Code
    
    Created by jnairn on Mon Nov 19 2001.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
*********************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "System/ArchiveData.hpp"
#include "Materials/MaterialBase.hpp"
#include "Custom_Tasks/CustomTask.hpp"
#include "Custom_Tasks/CalcJKTask.hpp"
#include "Custom_Tasks/PropagateTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Elements/ElementBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Exceptions/MPMTermination.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Exceptions/MPMWarnings.hpp"
#include "NairnMPM_Class/MPMTask.hpp"
#include "NairnMPM_Class/InitializationTask.hpp"
#include "NairnMPM_Class/MassAndMomentumTask.hpp"
#include "NairnMPM_Class/UpdateStrainsFirstTask.hpp"
#include "NairnMPM_Class/GridForcesTask.hpp"
#include "NairnMPM_Class/UpdateParticlesTask.hpp"
#include "NairnMPM_Class/UpdateStrainsLastTask.hpp"
#include "NairnMPM_Class/UpdateMomentaTask.hpp"
#include "NairnMPM_Class/RunCustomTasksTask.hpp"
#include "NairnMPM_Class/MoveCracksTask.hpp"
#include "NairnMPM_Class/ResetElementsTask.hpp"
#include <time.h>

// global analysis object
NairnMPM *fmobj=NULL;
MPMTask *firstMPMTask;

// global variables
double timestep=1.;			// time per MPM step (sec)
double strainTimestep;		// half time step when in US_AVG method
double mtime=0.;			// time for current step (sec)
double propTime=1e15;		// time interval between propagation calculations (sec)
int maxCrackFields=1;		// Maximum crack velocity fields at a node in a material (it is MAX_FIELDS_FOR_CRACKS if there are cracks)
int maxMaterialFields;		// Maximum velocity fields or number of independent materials

#pragma mark CONSTRUCTORS

// Constructor
NairnMPM::NairnMPM()
{
	version=8;						// main version
	subversion=1;					// subversion (must be < 10)
	buildnumber=0;					// build number
	mpmApproach=USAVG_METHOD;		// mpm method
	ptsPerElement=4;				// number of points per element (2D default, 3D changes it to 8)
	propagate[0]=propagate[1]=NO_PROPAGATION;						// default crack propagation type
	propagateDirection[0]=propagateDirection[1]=DEFAULT_DIRECTION;	// default crack propagation direction
	propagateMat[0]=propagateMat[1]=0;								// default is new crack with no traction law
	hasTractionCracks=FALSE;		// if any crack segment has a traction law material
	maxtime=1.;						// maximum time (sec)
	FractCellTime=.5;				// fraction cell crossed in 1 step at wave spd
	mstep=0;						// step number
	volumeExtrap=FALSE;				// set if need volume extrapolations
	warnParticleLeftGrid=1;			// abort when this many leave the grid
	multiMaterialMode=false;		// multi-material mode
	hasRigidContactParticles=false;	// rigid contact particles in multimaterial mode
	
	// initialize objects
	archiver=new ArchiveData();		// archiving object
}

#pragma mark METHODS

// start analysis
void NairnMPM::StartAnalysis(bool abort)
{
	// Active Transport Tasks
	if(DiffusionTask::active)
	{	transportTasks=new DiffusionTask();
	}
	if(ConductionTask::active)
	{	conduction=new ConductionTask();
		if(transportTasks)
			transportTasks->nextTask=conduction;
		else
			transportTasks=conduction;
	}
	
	// start results file
	StartResultsOutput();
	
	// Do MPM analysis
	MPMAnalysis(abort);
}

// Do the MPM analysis
void NairnMPM::MPMAnalysis(bool abort)
{
    char fline[100];
    CustomTask *nextTask;
    
    //---------------------------------------------------
	// initialize
    time(&startTime);
	startCPU=clock();
	
	//---------------------------------------------------
	// Do Preliminary MPM Calculations
	PreliminaryCalcs();

	// finish isothermal ramp (now that have time step)
	thermal.SetParameters(timestep);

	//---------------------------------------------------
	// Create custom tasks

	// if there are cracks, create J/K task and optionally a propagation task
	//		(it is essential for propagation task to be after the JK task)
	//		(insert this tasks before other custom tasks)
	if(firstCrack!=NULL)
	{	if(propagate[0] || archiver->WillArchiveJK(FALSE))
		{	nextTask=new CalcJKTask();
			if(propagate[0])
			{   nextTask=new PropagateTask();
				theJKTask->nextTask=nextTask;
			}
			nextTask->nextTask=theTasks;
			theTasks=theJKTask;
			ElementBase::AllocateNeighbors();
		}
	}

	// see if any need initializing
	if(theTasks!=NULL)
	{	PrintSection("SCHEDULED CUSTOM TASKS");
		nextTask=theTasks;
		while(nextTask!=NULL)
			nextTask=nextTask->Initialize();
		cout << endl;
	}

	// need particle volume? (assumes only transport tasks do)
	if(transportTasks) volumeExtrap=TRUE;
    
	//---------------------------------------------------
	// Create all the step tasks
	
	// TASK 0: INITIALIZATION
	MPMTask *lastMPMTask,*nextMPMTask;
	lastMPMTask=firstMPMTask=(MPMTask *)new InitializationTask("Initialization");
	
	// TASK 1: MASS MATRIX
	nextMPMTask=(MPMTask *)new MassAndMomentumTask("Mass matrix and momentum extrapolation");
	lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
	lastMPMTask=nextMPMTask;
	
	// TASK 2: UPDATE STRAINS FIRST AND USAVG
    if(mpmApproach==USF_METHOD || mpmApproach==USAVG_METHOD)
	{	nextMPMTask=(MPMTask *)new UpdateStrainsFirstTask("Update strains first");
		lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
		lastMPMTask=nextMPMTask;
	}
	
	// TASK 3: FORCES
	nextMPMTask=(MPMTask *)new GridForcesTask("Get grid forces");
	lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
	lastMPMTask=nextMPMTask;
    
	// TASK 4: UPDATE MOMENTA
	nextMPMTask=(MPMTask *)new UpdateMomentaTask("Update momenta");
	lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
	lastMPMTask=nextMPMTask;
    
	// TASK 5: UPDATE PARTICLES
	nextMPMTask=(MPMTask *)new UpdateParticlesTask("Update particles");
	lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
	lastMPMTask=nextMPMTask;
	
	// TASK 6: UPDATE STRAINS LAST AND USAVG
	if(mpmApproach==SZS_METHOD || mpmApproach==USAVG_METHOD)
	{	nextMPMTask=(MPMTask *)new UpdateStrainsLastTask("Update strains last");
		lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
		lastMPMTask=nextMPMTask;
	}
	
	// TASK 7: CUSTOM TASKS
	if(theTasks!=NULL)
	{	nextMPMTask=(MPMTask *)new RunCustomTasksTask("Run custom tasks");
		lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
		lastMPMTask=nextMPMTask;
	}
	
	// TASK 8a: MOVE CRACKS
	if(firstCrack!=NULL)
	{	nextMPMTask=(MPMTask *)new MoveCracksTask("Move cracks");
		lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
		lastMPMTask=nextMPMTask;
	}
	
	// TASK 8b: RESET ELEMEMTS
	nextMPMTask=(MPMTask *)new ResetElementsTask("Reset elements");
	lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
	lastMPMTask=nextMPMTask;
	
	try
	{	//---------------------------------------------------
		// Archiving
		archiver->BeginArchives(IsThreeD());
		archiver->ArchiveResults(mtime);
		
		// optional validation of parameters
		ValidateOptions();
		
		// exit if do not want analysis
		if(abort) mtime=maxtime+1;
		
		// ------------------------------------------------------
		// Main MPM Loop
		while(mtime<=maxtime)
		{	// MPM Calculations
			mstep++;			// step number
            MPMStep();
		
			// advance time and archive if desired
			mtime+=timestep;
			archiver->ArchiveResults(mtime);
        }
	}
	catch(MPMTermination term)
	{	// calculation stopped, but still report results
		mtime+=timestep;
		archiver->ForceArchiving();
		archiver->ArchiveResults(mtime);
		cout << endl;
		PrintSection("ABNORMAL TERMINATION");
		term.Display(mstep,mtime);
	}
	catch(CommonException err)
	{	// fatal error - exit to main
		cout << endl;
		PrintSection("ABNORMAL ERROR");
		throw err;
	}
	catch(const char *errMsg)
	{	// string error - exit to main
		throw errMsg;
	}
	catch(...)
	{	// unknown error - exit to main
		throw "Unknown exception in MPMStep() in MPM Loop in NairnMPM.cpp";
	}
    cout << endl;
	
    //---------------------------------------------------
    // Report warnings
	warnings.Report();

    //---------------------------------------------------
    // Report on time of exectution
    double execTime=ElapsedTime();						// elpased time in secs
	double cpuTime=CPUTime();							// cpu time in secs
	
    PrintSection("EXECUTION TIMES AND MEMORY");
    cout << "Calculation Steps: " << mstep << endl;
    
    sprintf(fline,"Elapsed Time: %.3lf secs\n",execTime);
    cout << fline;
    
    sprintf(fline,"CPU Time: %.3lf secs\n",cpuTime);
    cout << fline;
	
	if(mstep>0)
	{	sprintf(fline,"Elapsed Time per Step: %.3lf ms\n",1000.*execTime/((double)mstep));
		cout << fline;
		
		double timePerStep=1000.*cpuTime/((double)mstep);
		sprintf(fline,"CPU Time per Step: %.3lf ms\n",timePerStep);
		cout << fline;
		
#ifdef _PROFILE_TASKS_
		MPMTask *nextMPMTask=firstMPMTask;
		while(nextMPMTask!=NULL)
		{	nextMPMTask->WriteProfileResults(mstep,timePerStep);
			nextMPMTask=(MPMTask *)nextMPMTask->GetNextTask();
		}
#endif
	}
    
    //---------------------------------------------------
    // Trailer
    cout << "\n***** NairnMPM RUN COMPLETED\n";
}

// Main analysis loop for MPM analysis
// Made up of tasks created in MPMAnalysis()
void NairnMPM::MPMStep(void)
{
	// Step initialization
#ifdef LOG_PROGRESS
	char logLine[200];
	archiver->ClearLogFile();
	sprintf(logLine,"Step #%d: Initialize",mstep);
	archiver->WriteLogFile(logLine,NULL);
#endif
	
	// loop through the tasks
	MPMTask *nextMPMTask=firstMPMTask;
	while(nextMPMTask!=NULL)
	{
#ifdef LOG_PROGRESS
		nextMPMTask->WriteLogFile();
#endif
		nextMPMTask->Execute();
		nextMPMTask=(MPMTask *)nextMPMTask->GetNextTask();
#ifdef LOG_PROGRESS
		archiver->WriteLogFile("           Done",NULL);
#endif
	}

}

/**********************************************************
	Preliminary MPM Calclations prior to analysis
	1. Loop over particles
		- Check time step vs element size
		- check propagation time step
		- zero external forces on particle
		- Initialize history variables
	2. Loop over elements
	3. Print information about particles
**********************************************************/

void NairnMPM::PreliminaryCalcs(void)
{
    int p,i;
    short matid;
    double area,volume,rho,crot,tst,tmin=1e15;
    double dcell;
    char fline[200];
	
	// Loop over elements, if needed, to determine type of grid
	if(mpmgrid.GetCartesian()==UNKNOWN_GRID)
	{	int userCartesian=FALSE;
		double dx,dy,dz,gridx=0.,gridy=0.,gridz=0.;
		for(i=0;i<nelems;i++)
		{	if(!theElements[i]->Orthogonal(&dx,&dy,&dz))
			{	userCartesian=FALSE;
				break;
			}
			if(!userCartesian)
			{	// first element, set grid size (dz=0 in 2D) and set userCartesion flag
				gridx=dx;
				gridy=dy;
				gridz=dz;
				userCartesian=TRUE;
			}
			else
			{	// on sebsequent elements, if size does not match current values, than not Cartesian grod so give up
				if(!DbleEqual(gridx,dx) || !DbleEqual(gridy,dy) || !DbleEqual(gridz,dz))
				{	userCartesian=FALSE;
					break;
				}
			}
		}
		mpmgrid.SetCartesian(userCartesian,gridx,gridy,gridz);
	}
	
	// any material type checks needed?
	for(p=0;p<nmat;p++)
		theMaterials[p]->PreliminaryMatCalcs();
	
	// future - make this a parameter that can be input
	double PropFractCellTime=FractCellTime;
	double minSize=mpmgrid.GetMinCellDimension()/10.;	// in cm
    
    // loop over material points
	maxMaterialFields=0;
    for(p=0;p<nmpms;p++)
	{	matid=mpm[p]->MatID();
		if(matid>=nmat)
			throw CommonException("Material point with an undefined material type","NairnMPM::PreliminaryCalcs");
		if(theMaterials[matid]->isTractionLaw())
			throw CommonException("Material point with traction-law material","NairnMPM::PreliminaryCalcs");
		maxMaterialFields=theMaterials[matid]->SetField(maxMaterialFields,multiMaterialMode,matid);
		if(maxMaterialFields<0)
			throw CommonException("Rigid material with contact not allowed in single material mode MPM","NairnMPM::PreliminaryCalcs");
		if(theMaterials[matid]->RigidBC()) continue;
	
		// element and mp properties
		if(IsThreeD())
		{	volume=theElements[mpm[p]->ElemID()]->GetVolume()/1000.;	// in cm^3
			dcell = (minSize>0.) ? minSize : pow(volume,1./3.) ;
		}
		else
		{	area=theElements[mpm[p]->ElemID()]->GetArea()/100.;	// in cm^2
			volume=mpm[p]->thickness()*area/10.;				// in cm^2
			dcell = (minSize>0.) ? minSize : sqrt(area) ;
		}
		rho=theMaterials[matid]->rho;					// in g/cm^3
        
        // assumes same number of points for all elements
		mpm[p]->InitializeMass(rho*volume/((double)ptsPerElement));			// in g
		
		// done if rigid contact material in multimaterial mode (mass will be in mm^3 and will be particle volume)
		if(theMaterials[matid]->Rigid())
		{	hasRigidContactParticles=true;
			continue;
		}
        
        // check time step
        crot=theMaterials[matid]->WaveSpeed(IsThreeD())/10.;		// in cm/sec
		tst=FractCellTime*dcell/crot;					// in sec
        if(tst<tmin) tmin=tst;
        
        // propagation time (in sec)
        tst=PropFractCellTime*dcell/crot;
        if(tst<propTime) propTime=tst;
        
        // zero external forces on this particle
		ZeroVector(mpm[p]->GetPFext());
        
        // initialize history-dependent material data on this particle
        mpm[p]->SetHistoryPtr(theMaterials[matid]->MaterialData());
		
		// concentration potential
		if(mpm[p]->pConcentration<0.)
		{	double potential=-mpm[p]->pConcentration/theMaterials[matid]->concSaturation;
			if(potential>1.000001)
				throw CommonException("Material point with concentration potential > 1","NairnMPM::PreliminaryCalcs");
			if(potential>1.) potential=1.;
			mpm[p]->pConcentration=potential;
		}
		
		// Transport property time steps
		TransportTask *nextTransport=transportTasks;
		while(nextTransport!=NULL)
			nextTransport=nextTransport->TransportTimeStep(matid,dcell,&tmin);
	}
	
	// multimaterial checks and contact initialization
	if(maxMaterialFields==0)
		throw CommonException("No material points found with an actual material","NairnMPM::PreliminaryCalcs");
	else if(maxMaterialFields==1)
		multiMaterialMode=FALSE;
	else
		contact.MaterialContactPairs(maxMaterialFields);
	
    // verify time step and make smaller if needed
    if(tmin<timestep) timestep=tmin;
	strainTimestep = (mpmApproach==USAVG_METHOD) ? timestep/2. : timestep ;
	
	// propagation time step (no less than timestep)
    if(propTime<timestep) propTime=timestep;
	
    // Print particle information oand other preliminary calc results
    PrintSection("FULL MASS MATRIX");
    
    sprintf(fline,"Number of Material Points: %d",nmpms);
    cout << fline << endl;
	
	// background grid info
	mpmgrid.Output(ptsPerElement);
    
    sprintf(fline,"Adjusted time step (ms): %.7e",1000.*timestep);
    cout << fline << endl;
    
	// contact law materials and cracks
	contact.SetNormalCODCutoff(mpmgrid.GetMinCellDimension());
	
    // progation time step and other settings when has cracks
    if(firstCrack!=NULL)
	{	if(propagate[0])
		{   sprintf(fline,"Propagation time step (ms): %.7e",1000.*propTime);
			cout << fline << endl;
		}
		
		// let all cracks do prelimnary calcs
		CrackHeader *nextCrack=firstCrack;
		while(nextCrack!=NULL)
		{	nextCrack->PreliminaryCrackCalcs();
			if(nextCrack->GetHasTractionLaws())
				hasTractionCracks=TRUE;
			nextCrack=(CrackHeader *)nextCrack->GetNextObject();
		}
		
		maxCrackFields=MAX_FIELDS_FOR_CRACKS;
		
		// warnings
		CrackHeader::warnThreeFields=warnings.CreateWarning("crack node with three velocity fields",-1L);
		CrackHeader::warnNodeOnCrack=warnings.CreateWarning("mesh node on a crack",-1L);
		CrackHeader::warnThreeCracks=warnings.CreateWarning("node with three cracks or unexpected velocity fields",-1L);
	}
	
	// create warnings
	warnParticleLeftGrid=warnings.CreateWarning("particle left the grid",warnParticleLeftGrid);
	
	// nodal point calculations
	NodalPoint::PreliminaryCalcs();
    
    // blank line
    cout << endl;
}

// Called just before time steps start
// Can insert code here to black runs with invalid options
void NairnMPM::ValidateOptions(void)
{	
	if(ElementBase::useGimp)
	{	if(!mpmgrid.CanDoGIMP())
			throw CommonException("GIMP not allowed unless using a generated regular mesh","NairnMPM::ValidateOptions");
		if(ptsPerElement!=4 && !IsThreeD())
			throw CommonException("GIMP requires 4 particles per element for 2D","NairnMPM::ValidateOptions");
		if(ptsPerElement!=8 && IsThreeD())
			throw CommonException("GIMP requires 8 particles per element for 3D","NairnMPM::ValidateOptions");
	}
	
	if(contact.hasImperfectInterface)
	{	if(mpmgrid.GetCartesian()<=0)
			throw CommonException("Imperfect interfaces require a cartesian mesh","NairnMPM::ValidateOptions");
	}
	
	if(IsThreeD())
	{	if(mpmgrid.GetCartesian()!=CUBIC_GRID && mpmgrid.GetCartesian()!=ORTHOGONAL_GRID)
			throw CommonException("3D calculations require an orthogonal grid","NairnMPM::ValidateOptions");
		if(ptsPerElement!=1 && ptsPerElement!=8)
			throw CommonException("3D analysis requires 1 or 8 particles per cell","NairnMPM::ValidateOptions");
	}
	else
	{	if(ptsPerElement!=1 && ptsPerElement!=4)
			throw CommonException("2D analysis requires 1 or 4 particles per cell","NairnMPM::ValidateOptions");
	}
	
	if(multiMaterialMode)
	{	if(!mpmgrid.CanDoGIMP())
			throw CommonException("Multimaterial mode is not allowed unless using a generated regular mesh","NairnMPM::ValidateOptions");
	}
	
	// check each material type (if it is used)
	int i;
	for(i=0;i<nmat;i++)
	{	if(theMaterials[i]->GetField()>=0)
			theMaterials[i]->ValidateUse(np);
	}
}

#pragma mark ACCESSORS

// return name, caller should delete
char *NairnMPM::CodeName(void)
{
	char *name=new char[9];
	strcpy(name,"NairnMPM");
	return name;
}

// verify analysis type
bool NairnMPM::ValidAnalysisType(void)
{
	// change defaults for 3D - this will be called in <Header> which will be before any user changes
	if(np==THREED_MPM)
	{	ptsPerElement=8;				// number of points per element
		nfree=3;
	}
		
	return np>BEGIN_MPM_TYPES && np<END_MPM_TYPES;
}

// Explain usage of this program
void NairnMPM::Usage()
{
	CoutCodeVersion();
    cout << "    Expects Xerces 3.1.1 or newer" << endl;
    cout << "\nUsage:\n"
            "    NairnMPM [-options] <InputFile>\n\n"
            "This program reads the <InputFile> and does an MPM analysis.\n"
            "The results summary is directed to standard output. Errors are\n"
            "directed to standard error. The numerical results are saved\n"
            "to archived files as directed in <InputFile>.\n\n"
            "Options:\n"
            "    -a          Abort after setting up problem but before\n"
			"                   MPM steps begin. Initial conditions will\n"
			"                   be archived\n"
            "    -H          Show this help.\n"
            "    -r          Reverse byte order in archive files\n"
            "                   (default is to not reverse the bytes)\n"
            "    -v          Validate input file if DTD is provided in !DOCTYPE\n"
            "                   (default is to skip validation)\n\n"
            "See http://oregonstate.edu/nairnj for documentation.\n\n"
          <<  endl;
}

// elapsed CPU time (probably actual time)
double NairnMPM::ElapsedTime(void)
{	time_t currentTime;
    time(&currentTime);
    return (double)difftime(currentTime,startTime);
}

// elapsed CPU time (probably actual time)
double NairnMPM::CPUTime(void) { return (double)(clock()-startCPU)/CLOCKS_PER_SEC; }

// if crack develops tractionlaw, call here to turn it on
void NairnMPM::SetHasTractionCracks(bool setting) { hasTractionCracks=setting; }

