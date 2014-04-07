/*********************************************************************
    NairnMPM.cpp
    nairn-mpm-fea
    
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
#include "Boundary_Conditions/MatPtTractionBC.hpp"
#include "Boundary_Conditions/MatPtFluxBC.hpp"
#include "Boundary_Conditions/MatPtHeatFluxBC.hpp"
#include "Patches/GridPatch.hpp"
#include "Boundary_Conditions/NodalConcBC.hpp"
#include "Boundary_Conditions/NodalTempBC.hpp"
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
int maxMaterialFields;		// Maximum velocity fields or number of independent materials in multimaterial mode (=1 in single mat mode)
int numActiveMaterials;		// Number of non-rigid materials used by at least one material point
int maxShapeNodes=10;		// Maximum number of nodes for a particle (plus 1)

#pragma mark CONSTRUCTORS

// Constructor
NairnMPM::NairnMPM()
{
#ifdef _OSParticulas_
	version=1;						// main version
	subversion=0;					// subversion (must be < 10)
	buildnumber=0;					// build number
#else
	version=10;						// main version
	subversion=1;					// subversion (must be < 10)
	buildnumber=0;					// build number
#endif
	mpmApproach=USAVG_METHOD;		// mpm method
	ptsPerElement=4;				// number of points per element (2D default, 3D changes it to 8)
	propagate[0]=propagate[1]=NO_PROPAGATION;						// default crack propagation type
	propagateDirection[0]=propagateDirection[1]=DEFAULT_DIRECTION;	// default crack propagation direction
	propagateMat[0]=propagateMat[1]=0;								// default is new crack with no traction law
	hasTractionCracks=FALSE;		// if any crack segment has a traction law material
	maxtime=1.;						// maximum time (sec)
	FractCellTime=.5;				// fraction cell crossed in 1 step at wave speed (CFL convergence condition)
    PropFractCellTime=-1.;          // fracture cell crossed in 1 step for propagation time step (currently not user settable)
	mstep=0;						// step number
	warnParticleLeftGrid=-1;		// abort when this many leave the grid
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
	{	diffusion = new DiffusionTask();
		transportTasks = diffusion;
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
	// Do Preliminary MPM Calculations
	PreliminaryCalcs();

	//---------------------------------------------------
	// Create custom tasks

	// if there are cracks, create J/K task and optionally a propagation task
	//		(it is essential for propagation task to be after the JK task)
	//		(insert this task before other custom tasks)
	if(firstCrack!=NULL)
	{	if(propagate[0] || archiver->WillArchiveJK(FALSE))
		{   nextTask=new CalcJKTask();
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

	//---------------------------------------------------
	// Create all the step tasks
	
	// TASK 0: INITIALIZATION
	MPMTask *lastMPMTask,*nextMPMTask;
	lastMPMTask=firstMPMTask=(MPMTask *)new InitializationTask("Initialization");
	
	// TASK 1: MASS MATRIX
	nextMPMTask=(MPMTask *)new MassAndMomentumTask("Mass/momentum_extrapolation");
	lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
	lastMPMTask=nextMPMTask;
	
	// TASK 2: UPDATE STRAINS FIRST AND USAVG
    if(mpmApproach==USF_METHOD || mpmApproach==USAVG_METHOD)
	{	nextMPMTask=(MPMTask *)new UpdateStrainsFirstTask("Update_strains_first");
		lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
		lastMPMTask=nextMPMTask;
	}
	
	// TASK 3: FORCES
	nextMPMTask=(MPMTask *)new GridForcesTask("Grid_forces");
	lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
	lastMPMTask=nextMPMTask;
    
	// TASK 4: UPDATE MOMENTA
	nextMPMTask=(MPMTask *)new UpdateMomentaTask("Update_momenta");
	lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
	lastMPMTask=nextMPMTask;
    
	// TASK 5: UPDATE PARTICLES
	nextMPMTask=(MPMTask *)new UpdateParticlesTask("Update_particles");
	lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
	lastMPMTask=nextMPMTask;
	
	// TASK 6: UPDATE STRAINS LAST AND USAVG
	if(mpmApproach==SZS_METHOD || mpmApproach==USAVG_METHOD)
	{	nextMPMTask=(MPMTask *)new UpdateStrainsLastTask("Update_strains_last");
		lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
		lastMPMTask=nextMPMTask;
	}
	
	// TASK 7: CUSTOM TASKS
	if(theTasks!=NULL)
	{	nextMPMTask=(MPMTask *)new RunCustomTasksTask("Run_custom_tasks");
		lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
		lastMPMTask=nextMPMTask;
	}
	
	// TASK 8a: MOVE CRACKS
	if(firstCrack!=NULL)
	{	nextMPMTask=(MPMTask *)new MoveCracksTask("Move_cracks");
		lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
		lastMPMTask=nextMPMTask;
	}
	
	// TASK 8b: RESET ELEMEMTS
	nextMPMTask=(MPMTask *)new ResetElementsTask("Reset_elements");
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
	catch(CommonException term)
	{	// calculation stopped, but still report results
		mtime+=timestep;
		archiver->ForceArchiving();
		archiver->ArchiveResults(mtime);
		cout << endl;
		PrintSection("ABNORMAL TERMINATION");
		term.Display(mstep,mtime);
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
	{	double eTimePerStep = 1000.*execTime/((double)mstep);
		sprintf(fline,"Elapsed Time per Step: %.3lf ms\n",eTimePerStep);
		cout << fline;
		
		double timePerStep=1000.*cpuTime/((double)mstep);
		sprintf(fline,"CPU Time per Step: %.3lf ms\n",timePerStep);
		cout << fline;
		
#ifdef _PROFILE_TASKS_
		MPMTask *nextMPMTask=firstMPMTask;
		while(nextMPMTask!=NULL)
		{	nextMPMTask->WriteProfileResults(mstep,timePerStep,eTimePerStep);
			nextMPMTask=(MPMTask *)nextMPMTask->GetNextTask();
		}
#endif
	}
    
    //---------------------------------------------------
    // Trailer
    cout << "\n***** " << CodeName() << " RUN COMPLETED\n";
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
#ifdef _PROFILE_TASKS_
		double beginTime=fmobj->CPUTime();
		double beginETime=fmobj->ElapsedTime();
#endif
		nextMPMTask->Execute();
#ifdef _PROFILE_TASKS_
		nextMPMTask->TrackTimes(beginTime,beginETime);
#endif
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
	4. Initialize thermal calculations
**********************************************************/

void NairnMPM::PreliminaryCalcs(void)
{
    int p,i;
    short matid;
    double area,volume,rho,crot,tst,tmin=1e15;
    double dcell;
    char fline[200];
    
    // are both the system and the particles isolated?
    if(!ConductionTask::active && ConductionTask::IsSystemIsolated())
    {   MaterialBase::isolatedSystemAndParticles = TRUE;
    }
    
	// Loop over elements, if needed, to determine type of grid
	if(mpmgrid.GetCartesian()==UNKNOWN_GRID)
	{	int userCartesian=FALSE;
		double dx,dy,dz,gridx=0.,gridy=0.,gridz=0.;
		for(i=0;i<nelems;i++)
		{	if(!theElements[i]->Orthogonal(&dx,&dy,&dz))
            {   // exit if find one that is not orthongal
				userCartesian=FALSE;
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
			{	// on sebsequent elements, if size does not match current values, than not Cartesian grid so give up
				if(!DbleEqual(gridx,dx) || !DbleEqual(gridy,dy) || !DbleEqual(gridz,dz))
				{	userCartesian=FALSE;
					break;
				}
			}
		}
		mpmgrid.SetCartesian(userCartesian,gridx,gridy,gridz);
	}
    
    // CPDI factors if needed
    ElementBase::InitializeCPDI(IsThreeD());
	
    // future - make PropFractCellTime a user parameter, which not changed here if user picked it
    if(PropFractCellTime<0.) PropFractCellTime=FractCellTime;
	double minSize=mpmgrid.GetMinCellDimension()/10.;                   // in cm
    
    // loop over material points
	maxMaterialFields = 0;
	numActiveMaterials = 0;
    nmpmsNR = 0;
	int firstRigidPt = -1;
    for(p=0;p<nmpms;p++)
	{	// verify material is defined and sets if field number (in in multimaterial mode)
		matid=mpm[p]->MatID();
		if(matid>=nmat)
			throw CommonException("Material point with an undefined material type","NairnMPM::PreliminaryCalcs");
		if(theMaterials[matid]->isTractionLaw())
			throw CommonException("Material point with traction-law material","NairnMPM::PreliminaryCalcs");
		maxMaterialFields=theMaterials[matid]->SetField(maxMaterialFields,multiMaterialMode,matid,numActiveMaterials);
		
		// nothing left if rigid material is a BC rigid material
		if(theMaterials[matid]->RigidBC())
		{	if(firstRigidPt<0) firstRigidPt=p;
			continue;
		}
	
		// element and mp properties
		if(IsThreeD())
		{	volume=theElements[mpm[p]->ElemID()]->GetVolume()/1000.;	// in cm^3
			dcell = (minSize>0.) ? minSize : pow(volume,1./3.) ;
		}
		else
		{	// when axisymmetric, thickness is particle radial position, which gives mp = rho*Ap*Rp
			area=theElements[mpm[p]->ElemID()]->GetArea()/100.;	// in cm^2
			volume=mpm[p]->thickness()*area/10.;				// in cm^2
			dcell = (minSize>0.) ? minSize : sqrt(area) ;
		}
		rho=theMaterials[matid]->rho;					// in g/cm^3
        
        // assumes same number of points for all elements
        // for axisyymmeric xp = rho*Ap*volume/(# per element)
		mpm[p]->InitializeMass(rho*volume/((double)ptsPerElement));			// in g
		
		// done if rigid contact material in multimaterial mode
        // mass will be in mm^3 and will be particle volume
		if(theMaterials[matid]->Rigid())
		{	hasRigidContactParticles=true;
			if(firstRigidPt<0) firstRigidPt=p;
            
            // CPDI domain data
            if(!mpm[p]->AllocateCPDIStructures(ElementBase::useGimp,IsThreeD()))
                throw CommonException("Out of memory allocating CPDI domain structures","NairnMPM::PreliminaryCalcs");
			continue;
		}
        
        // now a nonrigid particle
        nmpmsNR = p+1;
        
        // check time step
        crot=theMaterials[matid]->WaveSpeed(IsThreeD(),mpm[p])/10.;		// in cm/sec
		tst=FractCellTime*dcell/crot;                                   // in sec
        if(tst<tmin) tmin=tst;
        
        // propagation time (in sec)
        tst=PropFractCellTime*dcell/crot;
        if(tst<propTime) propTime=tst;
        
        // zero external forces on this particle
		ZeroVector(mpm[p]->GetPFext());
        
        // initialize history-dependent material data on this particle
        mpm[p]->SetHistoryPtr(theMaterials[matid]->InitHistoryData());
		
		// concentration potential
		if(mpm[p]->pConcentration<0.)
		{	double potential=-mpm[p]->pConcentration/theMaterials[matid]->concSaturation;
			if(potential>1.000001)
				throw CommonException("Material point with concentration potential > 1","NairnMPM::PreliminaryCalcs");
			if(potential>1.) potential=1.;
			mpm[p]->pConcentration=potential;
		}
        
        // material dependent initialization
        theMaterials[matid]->SetInitialParticleState(mpm[p],np);
		
		// Transport property time steps
		int numTransport = 0;
		TransportTask *nextTransport=transportTasks;
		while(nextTransport!=NULL)
		{	numTransport++;
			nextTransport=nextTransport->TransportTimeStep(matid,dcell,&tmin);
		}
        
        // CPDI domain data
        if(!mpm[p]->AllocateCPDIStructures(ElementBase::useGimp,IsThreeD()))
            throw CommonException("Out of memory allocating CPDI domain structures","NairnMPM::PreliminaryCalcs");
		
	}
	
	// reorder NonRigid followed by (Rigid Contact and Rigid BCs intermized)
	if(firstRigidPt<nmpmsNR && firstRigidPt>=0)
	{	p = firstRigidPt;
		
		// loop until reach end of non-rigid materials
		while(p<nmpmsNR)
		{	matid=mpm[p]->MatID();
			if(theMaterials[matid]->Rigid())
			{	// if rigid particle, switch with particle at nmpmsNR-1
				MPMBase *temp = mpm[p];
				mpm[p] = mpm[nmpmsNR-1];		// move last nonrigid (#nmpmsNR) to position p (p+1)
				mpm[nmpmsNR-1] = temp;			// move rigid particle (p+1) to rigid domain (#nmpmNR)
				nmpmsNR--;						// now 0-based pointer of previous NR particle
				
				// fix particle based boundary conditions
				ReorderPtBCs(firstLoadedPt,p,nmpmsNR);
				ReorderPtBCs(firstTractionPt,p,nmpmsNR);
				ReorderPtBCs(firstFluxPt,p,nmpmsNR);
				ReorderPtBCs(firstHeatFluxPt,p,nmpmsNR);
				
				// back up to new last nonrigid
				while(nmpmsNR>=0)
				{	matid=mpm[nmpmsNR-1]->MatID();
					if(!theMaterials[matid]->Rigid()) break;
					nmpmsNR--;
				}
				
			}
			
			// next particle
			p++;
		}
	}
	
	// reorder rigid particles as (Rigid Contact, Rigid BCs)
	nmpmsRC = nmpmsNR;
	if(hasRigidContactParticles)
	{	// back up to new last rigid contact particle
		nmpmsRC = nmpms;
		while(nmpmsRC > nmpmsNR)
		{	matid = mpm[nmpmsRC-1]->MatID();
			if(!theMaterials[matid]->RigidBC()) break;
			nmpmsRC--;
		}
		
		// loop until reach end of rigid contact
		p = nmpmsNR;
		while(p < nmpmsRC)
		{	matid=mpm[p]->MatID();
			if(theMaterials[matid]->RigidBC())
			{	// if rigid BC particle, switch with rigid contact particle at nmpmsRC-1
				MPMBase *temp = mpm[p];
				mpm[p] = mpm[nmpmsRC-1];
				mpm[nmpmsRC-1] = temp;
				nmpmsRC--;
				
				// back up to new last nonrigid
				while(nmpmsRC > nmpmsNR)
				{	matid=mpm[nmpmsRC-1]->MatID();
					if(!theMaterials[matid]->RigidBC()) break;
					nmpmsRC--;
				}
			}
			
			// next particle
			p++;
		}
	}
    
	// multimaterial checks and contact initialization
	if(maxMaterialFields==0 || numActiveMaterials==0)
		throw CommonException("No material points found with an actual material","NairnMPM::PreliminaryCalcs");
	else if(maxMaterialFields==1)
		multiMaterialMode=FALSE;
	else
		contact.MaterialContactPairs(maxMaterialFields);
	
    // ghost particles will need to now this setting when creating patches
    if(firstCrack!=NULL) maxCrackFields=MAX_FIELDS_FOR_CRACKS;
    
    // verify time step and make smaller if needed
    if(tmin<timestep) timestep=tmin;
	strainTimestep = (mpmApproach==USAVG_METHOD) ? timestep/2. : timestep ;
	
	// propagation time step (no less than timestep)
    if(propTime<timestep) propTime=timestep;
	
	// create patches or a single patch
	patches = mpmgrid.CreatePatches(np,numProcs);
    if(patches==NULL)
		throw CommonException("Out of memory creating the patches","NairnMPM::PreliminaryCalcs");
    
    // create buffers for copies of material properties
    UpdateStrainsFirstTask::CreatePropertyBuffers(GetTotalNumberOfPatches());
	
	//---------------------------------------------------
	// Finish Results File
	BoundaryCondition *nextBC;
	
    //---------------------------------------------------
    // Loaded Material Points
    if(firstLoadedPt!=NULL)
    {   PrintSection("MATERIAL POINTS WITH EXTERNAL FORCES");
        cout << "Point   DOF ID     Load (N)     Arg (ms/ms^-1)  Function\n"
        << "----------------------------------------------------------\n";
        nextBC=(BoundaryCondition *)firstLoadedPt;
        while(nextBC!=NULL)
            nextBC=nextBC->PrintBC(cout);
        cout << endl;
    }
	
	//---------------------------------------------------
    // Traction Loaded Material Points
    if(firstTractionPt!=NULL)
    {   PrintSection("MATERIAL POINTS WITH TRACTIONS");
        cout << "Point   DOF Face ID   Stress (MPa)    Arg (ms/ms^-1)  Function\n"
        << "----------------------------------------------------------------\n";
        nextBC=(BoundaryCondition *)firstTractionPt;
        while(nextBC!=NULL)
            nextBC=nextBC->PrintBC(cout);
        cout << endl;
    }
	
	//---------------------------------------------------
    // Diffusion boundary conditions
	if(DiffusionTask::active)
	{   PrintSection("NODAL POINTS WITH FIXED CONCENTRATIONS");
		cout << "  Node  ID    Conc (/csat)   Arg (ms/ms^-1)  Function\n"
		<< "------------------------------------------------------\n";
		nextBC=firstConcBC;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
		
		//---------------------------------------------------
		// Concentration Flux Material Points
		PrintSection("MATERIAL POINTS WITH CONCENTRATION FLUX");
		cout << " Point  DOF Face ID   Flux (mm/sec)   Arg (ms/ms^-1)  Function\n"
		<< "---------------------------------------------------------------\n";
		nextBC=(BoundaryCondition *)firstFluxPt;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
	}
	
	//---------------------------------------------------
    // Conduction boundary conditions
	if(ConductionTask::active)
	{   PrintSection("NODAL POINTS WITH FIXED TEMPERATURES");
		cout << " Node   ID   Temp (-----)   Arg (ms/ms^-1)  Function\n"
		<< "------------------------------------------------------\n";
		nextBC=firstTempBC;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
		
		//---------------------------------------------------
		// Heat Flux Material Points
		PrintSection("MATERIAL POINTS WITH HEAT FLUX");
		cout << " Point  DOF Face ID   Flux (W/m^2)    Arg (ms/ms^-1)  Function\n"
		<< "---------------------------------------------------------------\n";
		nextBC=(BoundaryCondition *)firstHeatFluxPt;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
	}

    // Print particle information oand other preliminary calc results
    PrintSection("FULL MASS MATRIX");
    
    sprintf(fline,"Number of Material Points: %d",nmpms);
    cout << fline << endl;
	
	// background grid info
	mpmgrid.Output(ptsPerElement,IsAxisymmetric());
    
    sprintf(fline,"Adjusted time step (ms): %.7e",1000.*timestep);
    cout << fline << endl;
    
    // prpagation time step and other settings when has cracks
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
		
		// warnings
		CrackHeader::warnNodeOnCrack=warnings.CreateWarning("mesh node on a crack",-1,5);
		CrackHeader::warnThreeCracks=warnings.CreateWarning("node with three cracks or unexpected velocity fields",-1,0);
	}
	
	// create warnings
	if(warnParticleLeftGrid<0)
	{	// use default setting to quit if 1% of particle leave the grid
		warnParticleLeftGrid = nmpms/100;
	}
	warnParticleLeftGrid=warnings.CreateWarning("particle has left the grid warning",warnParticleLeftGrid,0);
	
	// nodal point calculations
	NodalPoint::PreliminaryCalcs();
	
	// finish isothermal ramp (now that have time step)
	thermal.SetParameters(timestep);
    
    // blank line
    cout << endl;
}

// When NR particle p2 moves to p1, reset any point-based BCs that use that point
void NairnMPM::ReorderPtBCs(MatPtLoadBC *firstBC,int p1,int p2)
{
	while(firstBC!=NULL)
		firstBC = firstBC->ReorderPtNum(p1,p2);
}

// Called just before time steps start
// Can insert code here to black runs with invalid options
void NairnMPM::ValidateOptions(void)
{	
    // GIMP and CPDI require regular
    //  and qCPDI not allowed in 3D
	if(ElementBase::useGimp != POINT_GIMP)
    {   // using a GIMP method
		if(!mpmgrid.CanDoGIMP())
			throw CommonException("GIMP not allowed unless using a generated regular mesh","NairnMPM::ValidateOptions");
        if(ElementBase::useGimp == QUADRATIC_CPDI)
        {   if(IsThreeD())
                throw CommonException("3D does not allow qCPDI shape functions; use lCPDI instead","NairnMPM::ValidateOptions");
        }
	}
    else
    {   // in Classic MPM or POINT_GIMP, cannot use traction BCs
        if(firstTractionPt!=NULL || firstFluxPt!=NULL || firstHeatFluxPt!=NULL)
			throw CommonException("Traction and flux boundary conditions require use of a GIMP MPM method.","NairnMPM::ValidateOptions");
    }
    
    // Imperfect interface requires cartensian grid
	if(contact.hasImperfectInterface)
	{	if(mpmgrid.GetCartesian()<=0)
			throw CommonException("Imperfect interfaces require a cartesian mesh","NairnMPM::ValidateOptions");
	}
	
    // 3D requires orthogonal grid and 1 or 8 particles per element
    // 2D requires 1 or 4 particles per element
	if(IsThreeD())
	{	if(mpmgrid.GetCartesian()!=CUBIC_GRID && mpmgrid.GetCartesian()!=ORTHOGONAL_GRID)
			throw CommonException("3D calculations require an orthogonal grid","NairnMPM::ValidateOptions");
		if(ptsPerElement!=1 && ptsPerElement!=8 && ptsPerElement!=27)
			throw CommonException("3D analysis requires 1 or 8 or 27 particles per cell","NairnMPM::ValidateOptions");
	}
	else
	{	if(ptsPerElement!=1 && ptsPerElement!=4 && ptsPerElement!=9 && ptsPerElement!=16 && ptsPerElement!=25)
			throw CommonException("2D analysis requires 1, 4, 9, 16, or 25 particles per cell","NairnMPM::ValidateOptions");
	}

    // Axisymmetric requirements and adjustments
    if(IsAxisymmetric())
    {   if(ElementBase::useGimp == POINT_GIMP)
        {   // require cartesian grid
            if(mpmgrid.GetCartesian()<=0)
            {   throw CommonException("Axisymmetric with Classic MPM requires anorthogonal grid","NairnMPM::ValidateOptions");
            }
        }
        else if(ElementBase::useGimp == UNIFORM_GIMP)
		{	ElementBase::useGimp = UNIFORM_GIMP_AS;
			ElementBase::analysisGimp = UNIFORM_GIMP_AS;
        }
        else if(ElementBase::useGimp == LINEAR_CPDI)
		{	ElementBase::useGimp = LINEAR_CPDI_AS;
			ElementBase::analysisGimp = LINEAR_CPDI_AS;
        }
		else
		{   throw CommonException("Axisymmetric does not allow qCPDI shape functions","NairnMPM::ValidateOptions");
		}
    }
	
    // Multimaterial mode requires a regular grid
	if(multiMaterialMode)
	{	if(!mpmgrid.CanDoGIMP())
			throw CommonException("Multimaterial mode is not allowed unless using a generated regular mesh","NairnMPM::ValidateOptions");
	}

#ifdef _OPENMP
	if(!mpmgrid.CanDoGIMP())
		throw CommonException("Multicore parallel code requires a generated regular mesh","NairnMPM::ValidateOptions");
#endif
	
	// check each material type (but only if it is used it at least one material point)
	int i;
	for(i=0;i<nmat;i++)
	{	if(theMaterials[i]->GetField()>=0)
			theMaterials[i]->ValidateForUse(np);
	}
}

#pragma mark ACCESSORS

// return name
const char *NairnMPM::CodeName(void) const
{
#ifdef _OSParticulas_
    return "OSParticulas";
#else
	return "NairnMPM";
#endif
}

// verify analysis type
bool NairnMPM::ValidAnalysisType(void)
{
	// change defaults for 3D - this will be called in <Header> which will be before any user changes
	if(np==THREED_MPM)
	{	ptsPerElement=8;				// number of points per element
		nfree=3;
		maxShapeNodes=28;				// increase if CPDI is used
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

// if crack develops tractionlaw, call here to turn it on
void NairnMPM::SetHasTractionCracks(bool setting) { hasTractionCracks=setting; }

// Get Courant-Friedrichs-Levy condition factor for convergence
void NairnMPM::SetCFLCondition(double factor) { FractCellTime = factor; }
double NairnMPM::GetCFLCondition(void) { return FractCellTime; }
double *NairnMPM::GetCFLPtr(void) { return &FractCellTime; }

// Get Courant-Friedrichs-Levy condition factor for convergence for propagation calculations
double NairnMPM::GetPropagationCFLCondition(void) { return PropFractCellTime; }


