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
#include "NairnMPM_Class/InitVelocityFieldsTask.hpp"
#include "NairnMPM_Class/ProjectRigidBCsTask.hpp"
#include "NairnMPM_Class/ExtrapolateRigidBCsTask.hpp"
#include "NairnMPM_Class/PostExtrapolationTask.hpp"
#include "NairnMPM_Class/SetRigidContactVelTask.hpp"
#include "NairnMPM_Class/MassAndMomentumTask.hpp"
#include "NairnMPM_Class/UpdateStrainsFirstTask.hpp"
#include "NairnMPM_Class/GridForcesTask.hpp"
#include "NairnMPM_Class/PostForcesTask.hpp"
#include "NairnMPM_Class/UpdateParticlesTask.hpp"
#include "NairnMPM_Class/UpdateStrainsLastTask.hpp"
#include "NairnMPM_Class/UpdateStrainsLastContactTask.hpp"
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
#include "System/UnitsController.hpp"
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
// throws std::bad_alloc
NairnMPM::NairnMPM()
{
	version=11;						// main version
	subversion=4;					// subversion (must be < 10)
	buildnumber=0;					// build number

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
	skipPostExtrapolation=false;	// optionally do not extrapolate for post update strain updates
	plusParticleSpin=false;			// add particle spin feature
	
	// initialize objects
	archiver=new ArchiveData();		// archiving object
}

#pragma mark METHODS

// start analysis
// throws std::bad_alloc
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
	else
	{	// these can only be on when conduction is active
		ConductionTask::crackContactHeating = false;
		ConductionTask::matContactHeating = false;
		ConductionTask::crackTipHeating = false;
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
    
	//---------------------------------------------------
	// Do Preliminary MPM Calculations
	PreliminaryCalcs();

	//---------------------------------------------------
	// Create step tasks
	CreateTasks();
	
	try
	{	//---------------------------------------------------
		// Archiving
		if(!archiver->BeginArchives(IsThreeD(),maxMaterialFields))
			throw "No archiving was specified or multiple archiving blocks not monotonically increasing in start time";
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
	catch(CommonException& term)
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
		
        // profile task results
		MPMTask *nextMPMTask=firstMPMTask;
		while(nextMPMTask!=NULL)
		{	nextMPMTask->WriteProfileResults(mstep,timePerStep,eTimePerStep);
			nextMPMTask=(MPMTask *)nextMPMTask->GetNextTask();
		}
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
		double beginTime=fmobj->CPUTime();
		double beginETime=fmobj->ElapsedTime();
		nextMPMTask->Execute();
		nextMPMTask->TrackTimes(beginTime,beginETime);
        
        // on to next task
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
 
	throws CommonException()
**********************************************************/

void NairnMPM::PreliminaryCalcs(void)
{
    int p,i;
    short matid;
    double area,volume,rho,crot,tst,tmin=1e15;
    char fline[200];
	
    // are both the system and the particles isolated?
    if(!ConductionTask::active && ConductionTask::IsSystemIsolated())
    {   MaterialBase::isolatedSystemAndParticles = TRUE;
    }
    
	// only allows grids created with the Grid command
	// (note: next two section never occur uless support turned back on)
	if(mpmgrid.GetCartesian()==UNKNOWN_GRID)
		throw CommonException("Support for iunstructured grids is currently not available","NairnMPM::PreliminaryCalcs");
	
	// Loop over elements, if needed, to determine type of grid
	if(mpmgrid.GetCartesian()==UNKNOWN_GRID)
	{	int userCartesian = NOT_CARTESIAN;
		double dx,dy,dz,gridx=0.,gridy=0.,gridz=0.;
		for(i=0;i<nelems;i++)
		{	if(!theElements[i]->Orthogonal(&dx,&dy,&dz))
            {   // exit if find one that is not orthogonal
				userCartesian = NOT_CARTESIAN;
				break;
			}
			if(userCartesian == NOT_CARTESIAN)
			{	// first element, set grid size (dz=0 in 2D) and set userCartesion flag
				gridx=dx;
				gridy=dy;
				gridz=dz;
				userCartesian = SQUARE_GRID;
			}
			else if(userCartesian==SQUARE_GRID)
			{	// on sebsequent elements, if size does not match current values, than elements differ
                // in size. Set to variable grid but keep going to check if rest are orthogonal
				if(!DbleEqual(gridx,dx) || !DbleEqual(gridy,dy) || !DbleEqual(gridz,dz))
				{	if(IsThreeD())
					{   gridz = 1.;
						userCartesian = VARIABLE_ORTHOGONAL_GRID;
					}
					else
						userCartesian = VARIABLE_RECTANGULAR_GRID;
                }
			}
		}
		
		// unstructer grid - set type, but elements not set so meshInfo will not know many things
		// and horiz will lbe <=0 
		mpmgrid.SetCartesian(userCartesian,gridx,gridy,gridz);
	}
	
	// only allows orthogonal grids
	if(mpmgrid.GetCartesian()<=0 || mpmgrid.GetCartesian()==NOT_CARTESIAN_3D)
		throw CommonException("Support for non-cartesian grids is currently not available","NairnMPM::PreliminaryCalcs");

    // CPDI factors if needed
    ElementBase::InitializeCPDI(IsThreeD());
	
    // future - make PropFractCellTime a user parameter, which not changed here if user picked it
    if(PropFractCellTime<0.) PropFractCellTime=FractCellTime;
	double dcell = mpmgrid.GetMinCellDimension();                   // in mm
    
    // loop over material points
	maxMaterialFields = 0;
	numActiveMaterials = 0;
    nmpmsNR = 0;
	bool hasRigidContactParticles = false;			// this become local variable in this mode, when off it is class variable
	int firstRigidPt = -1;
    for(p=0;p<nmpms;p++)
	{	// verify material is defined and set its field number (if in multimaterial mode)
		matid=mpm[p]->MatID();
		if(matid>=nmat)
			throw CommonException("Material point with an undefined material type","NairnMPM::PreliminaryCalcs");
		
		// material point can't use traction law or contact law
		if(theMaterials[matid]->MaterialStyle()==TRACTION_MAT)
			throw CommonException("Material point with traction-law material","NairnMPM::PreliminaryCalcs");
		if(theMaterials[matid]->MaterialStyle()==CONTACT_MAT)
			throw CommonException("Material point with contact-law material","NairnMPM::PreliminaryCalcs");
		
        // initialize history-dependent material data on this particle
		// might need initialized history data so set it now, but nonrigid only
		if(!theMaterials[matid]->Rigid())
			mpm[p]->SetHistoryPtr(theMaterials[matid]->InitHistoryData(NULL,mpm[p]));
		
		// Set material field (always [0] and returns 1 in single material mode)
		// Also verify allowsCracks - can only be true if multimaterial mode; if rigid must be rigid contact material
		maxMaterialFields=theMaterials[matid]->SetField(maxMaterialFields,multiMaterialMode,matid,numActiveMaterials);
		
		// element and mp properties
		if(IsThreeD())
		{	volume=theElements[mpm[p]->ElemID()]->GetVolume();	// in mm^3
		}
		else
		{	// when axisymmetric, thickness is particle radial position, which gives mp = rho*Ap*Rp
			area=theElements[mpm[p]->ElemID()]->GetArea();		// in mm^2
			volume=mpm[p]->thickness()*area;					// in mm^2
		}
		rho=theMaterials[matid]->GetRho(mpm[p]);							// in g/mm^3
        
        // assumes same number of points for all elements (but subclass can override)
        // for axisymmetric xp = rho*Ap*volume/(# per element)
		mpm[p]->InitializeMass(rho,volume/((double)ptsPerElement),plusParticleSpin);			// in g
		
		// done if rigid - mass will be in mm^3 and will be particle volume
		if(theMaterials[matid]->Rigid())
		{	// Trap rigid BC (used to be above element and mp properties,
			//     but moved here to get volume)
			if(theMaterials[matid]->RigidBC())
			{	if(firstRigidPt<0) firstRigidPt=p;
				continue;
			}
			
			// rigid contact material in multimaterial mode
			// Note: multimaterial mode implied here because otherwise fatal error when SetField() above
			hasRigidContactParticles=true;
			if(firstRigidPt<0) firstRigidPt=p;

			// CPDI or GIMP domain data only for rigid contact particles (non-rigid done below)
			if(!mpm[p]->AllocateCPDIorGIMPStructures(ElementBase::useGimp,IsThreeD()))
				throw CommonException("Out of memory allocating CPDI domain structures","NairnMPM::PreliminaryCalcs");
				
			continue;
		}
        
        // now a nonrigid particle
        nmpmsNR = p+1;
        
        // zero external forces on this particle
		ZeroVector(mpm[p]->GetPFext());
        
		// concentration potential
		if(mpm[p]->pConcentration<0.)
		{	double potential=-mpm[p]->pConcentration/mpm[p]->GetConcSaturation();
			if(potential>1.000001)
				throw CommonException("Material point with concentration potential > 1","NairnMPM::PreliminaryCalcs");
			if(potential>1.) potential=1.;
			mpm[p]->pConcentration=potential;
		}
        
        // material dependent initialization
        theMaterials[matid]->SetInitialParticleState(mpm[p],np,0);
		
        // check time step
        crot=theMaterials[matid]->WaveSpeed(IsThreeD(),mpm[p]);			// in mm/sec
		tst=FractCellTime*dcell/crot;                                   // in sec
        if(tst<tmin) tmin=tst;
        
        // propagation time (in sec)
        tst=PropFractCellTime*dcell/crot;
        if(tst<propTime) propTime=tst;
        
		// Transport property time steps
		int numTransport = 0;
		TransportTask *nextTransport=transportTasks;
		while(nextTransport!=NULL)
		{	numTransport++;
			nextTransport=nextTransport->TransportTimeStep(matid,dcell,&tmin);
		}
        
        // CPDI orGIMP domain data for nonrigid particles
        if(!mpm[p]->AllocateCPDIorGIMPStructures(ElementBase::useGimp,IsThreeD()))
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
    if(firstCrack!=NULL)
    {   if(firstCrack->GetNextObject()!=NULL)
            maxCrackFields=MAX_FIELDS_FOR_CRACKS;
        else
            maxCrackFields=MAX_FIELDS_FOR_ONE_CRACK;
    }
    
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
        cout << "Point   DOF ID     Load (" << UnitsController::Label(FEAFORCE_UNITS) << ")     Arg ("
			<< UnitsController::Label(BCARG_UNITS) << ")  Function\n"
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
        cout << "Point   DOF Face ID   Stress (" << UnitsController::Label(PRESSURE_UNITS) << ")    Arg ("
			<< UnitsController::Label(BCARG_UNITS) << ")  Function\n"
			<< "----------------------------------------------------------------\n";
        nextBC=(BoundaryCondition *)firstTractionPt;
        while(nextBC!=NULL)
            nextBC=nextBC->PrintBC(cout);
        cout << endl;
    }
	
	//---------------------------------------------------
    // Diffusion boundary conditions
	if(DiffusionTask::active && firstConcBC!=NULL)
	{   PrintSection("NODAL POINTS WITH FIXED CONCENTRATIONS");
		cout << "  Node  ID    Conc (/csat)   Arg (" << UnitsController::Label(BCARG_UNITS) << ")  Function\n"
			<< "------------------------------------------------------\n";
		nextBC=firstConcBC;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
	}
	
	//---------------------------------------------------
	// Concentration Flux Material Points
	if(DiffusionTask::active && firstFluxPt!=NULL)
	{	PrintSection("MATERIAL POINTS WITH CONCENTRATION FLUX");
		cout << " Point  DOF Face ID Flux (" << UnitsController::Label(BCCONCFLUX_UNITS) << ") Arg ("
			<< UnitsController::Label(BCARG_UNITS) << ")  Function\n"
			<< "---------------------------------------------------------------\n";
		nextBC=(BoundaryCondition *)firstFluxPt;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
	}
	
	//---------------------------------------------------
    // Conduction boundary conditions
	if(ConductionTask::active && firstTempBC!=NULL)
	{   PrintSection("NODAL POINTS WITH FIXED TEMPERATURES");
		cout << " Node   ID   Temp (-----)   Arg (" << UnitsController::Label(BCARG_UNITS) << ")  Function\n"
			<< "------------------------------------------------------\n";
		nextBC=firstTempBC;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
	}
		
	//---------------------------------------------------
	// Heat Flux Material Points
	if(ConductionTask::active && firstHeatFluxPt!=NULL)
	{	PrintSection("MATERIAL POINTS WITH HEAT FLUX");
		cout << " Point  DOF Face ID   Flux (" << UnitsController::Label(BCHEATFLUX_UNITS) << ")    Arg ("
			<< UnitsController::Label(BCARG_UNITS) << ")  Function\n"
			<< "---------------------------------------------------------------\n";
		nextBC=(BoundaryCondition *)firstHeatFluxPt;
		while(nextBC!=NULL)
			nextBC=nextBC->PrintBC(cout);
		cout << endl;
	}
	
	// non-standard particle sizes
	archiver->ArchivePointDimensions();

    // Print particle information and other preliminary calc results
    PrintSection("FULL MASS MATRIX");
    
    sprintf(fline,"Number of Material Points: %d",nmpms);
    cout << fline << endl;
	
	// background grid info
	mpmgrid.Output(ptsPerElement,IsAxisymmetric());
    
    sprintf(fline,"Adjusted time step (%s): %.7e",UnitsController::Label(ALTTIME_UNITS),timestep*UnitsController::Scaling(1.e3));
    cout << fline << endl;
    
    // prpagation time step and other settings when has cracks
    if(firstCrack!=NULL)
	{	if(propagate[0])
		{   sprintf(fline,"Propagation time step (%s): %.7e",UnitsController::Label(ALTTIME_UNITS),propTime*UnitsController::Scaling(1.e3));
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
		CrackHeader::warnNodeOnCrack=warnings.CreateWarning("unexpected crack side; possibly caused by node or particle on a crack",-1,5);
		CrackHeader::warnThreeCracks=warnings.CreateWarning("node with three or more cracks",-1,5);
	}
	
	// create warnings
	if(warnParticleLeftGrid<0)
	{	// use default setting to quit if 1% of particle leave the grid
		warnParticleLeftGrid = nmpms/100;
	}
	warnParticleLeftGrid=warnings.CreateWarning("particle has left the grid warning",warnParticleLeftGrid,0);
	
	// nodal point calculations
	NodalPoint::PreliminaryCalcs();
	
	// finish warnings
	warnings.CreateWarningList();
	
    // blank line
    cout << endl;
}

// create all the tasks needed for current simulation
// throws std::bad_alloc
void NairnMPM::CreateTasks(void)
{
    CustomTask *nextTask;
	
	// if there are cracks, create J/K task and optionally a propagation task
	//		(it is essential for propagation task to be after the JK task)
	//		(insert these tasks before other custom tasks)
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
		else
		{   // Turn off any J-K settings that can cause problems
			JTerms = -1;
		}
	}
    else
    {   // Turn off any crack settings that can cause problems
        JTerms = -1;
    }
	
	// see if any custom or transport tasks need initializing
	if(theTasks!=NULL || transportTasks!=NULL)
	{	PrintSection("SCHEDULED CUSTOM TASKS");
		
		// transport tasks
		TransportTask *nextTransport = transportTasks;
		while(nextTransport!=NULL)
			nextTransport = nextTransport->Initialize();
		
		// custom tasks
		nextTask=theTasks;
		while(nextTask!=NULL)
			nextTask=nextTask->Initialize();
		
		cout << endl;
	}
	
	//---------------------------------------------------
	// Create all the step tasks
	
	// INITIALIZATION Tasks
	MPMTask *lastMPMTask,*nextMPMTask;
	lastMPMTask=firstMPMTask=(MPMTask *)new InitializationTask("Initialize");
	
	if(firstCrack!=NULL || maxMaterialFields>1)
	{	nextMPMTask=(MPMTask *)new InitVelocityFieldsTask("Decipher Crack and Material Fields");
		lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
		lastMPMTask=nextMPMTask;
	}
	
	// MASS AND MOMENTUM TASKS
	
	// if rigid BCs by extrapolation, extrapolate now
	if(nmpms>nmpmsRC && MaterialBase::extrapolateRigidBCs)
	{	nextMPMTask=(MPMTask *)new ExtrapolateRigidBCsTask("Rigid BCs by Extrapolation");
		lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
		lastMPMTask=nextMPMTask;
	}
	
	// if rigid contact, add task to set their velocities
	if(nmpmsRC > nmpmsNR)
	{	nextMPMTask=(MPMTask *)new SetRigidContactVelTask("Set Rigid Contact Velocities");
		lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
		lastMPMTask=nextMPMTask;
	}
	
	nextMPMTask=(MPMTask *)new MassAndMomentumTask("Extrapolate Mass and Momentum");
	lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
	lastMPMTask=nextMPMTask;
	
	// if rigid BCs by extrapolation, project to BCs
	if(nmpms>nmpmsRC && !MaterialBase::extrapolateRigidBCs)
	{	nextMPMTask=(MPMTask *)new ProjectRigidBCsTask("Rigid BCs by Projection");
		lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
		lastMPMTask=nextMPMTask;
	}
	
	nextMPMTask=(MPMTask *)new PostExtrapolationTask("Post Extrapolation Tasks");
	lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
	lastMPMTask=nextMPMTask;
	
	// UPDATE STRAINS FIRST AND USAVG
    if(mpmApproach==USF_METHOD || mpmApproach==USAVG_METHOD)
	{	nextMPMTask=(MPMTask *)new UpdateStrainsFirstTask("Update Strains First");
		lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
		lastMPMTask=nextMPMTask;
	}
	
	// EXTRAPOLATE FORCES
	nextMPMTask=(MPMTask *)new GridForcesTask("Extrapolate Grid Forces");
	lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
	lastMPMTask=nextMPMTask;
    
	nextMPMTask=(MPMTask *)new PostForcesTask("Post Force Extrapolation Tasks");
	lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
	lastMPMTask=nextMPMTask;
    
	// UPDATE MOMENTA
	nextMPMTask=(MPMTask *)new UpdateMomentaTask("Update Momenta");
	lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
	lastMPMTask=nextMPMTask;
    
	// UPDATE PARTICLES
	nextMPMTask=(MPMTask *)new UpdateParticlesTask("Update Particles");
	lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
	lastMPMTask=nextMPMTask;
	
	// UPDATE STRAINS LAST AND USAVG
	// Energy calcs suggest the contact method, which re-extrapolates, should always be used
	if(mpmApproach==SZS_METHOD || mpmApproach==USAVG_METHOD)
	{	if(fmobj->skipPostExtrapolation)
		{	nextMPMTask=(MPMTask *)new UpdateStrainsLastTask("Update Strains Last");
			lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
			lastMPMTask=nextMPMTask;
		}
		else
		{	nextMPMTask=(MPMTask *)new UpdateStrainsLastContactTask("Update Strains Last with Extrapolation");
			lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
			lastMPMTask=nextMPMTask;
		}
	}
	
	// CUSTOM TASKS (including J Integral and crack propagation)
	if(theTasks!=NULL)
    {   nextMPMTask=(MPMTask *)new RunCustomTasksTask("Run Custom Tasks");
		lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
		lastMPMTask=nextMPMTask;
	}
	
	// MOVE CRACKS
	if(firstCrack!=NULL)
	{	nextMPMTask=(MPMTask *)new MoveCracksTask("Move Cracks");
		lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
		lastMPMTask=nextMPMTask;
	}
	
	// RESET ELEMEMTS
	nextMPMTask=(MPMTask *)new ResetElementsTask("Reset Elements");
	lastMPMTask->SetNextTask((CommonTask *)nextMPMTask);
	lastMPMTask=nextMPMTask;

}

// When NR particle p2 moves to p1, reset any point-based BCs that use that point
void NairnMPM::ReorderPtBCs(MatPtLoadBC *firstBC,int p1,int p2)
{
	while(firstBC!=NULL)
		firstBC = firstBC->ReorderPtNum(p1,p2);
}

// Called just before time steps start
// Can insert code here to black runs with invalid options
// throws CommonException()
void NairnMPM::ValidateOptions(void)
{
	// Disable non-structured or variable element grid sizes
	if(!mpmgrid.IsStructuredEqualElementsGrid())
	{	throw CommonException("Non-structured grids or grids with unequal element sizes are currently disabled in NairnMPM because they are not verified for all features.",
							  "NairnMPM::ValidateOptions");
	}
	
    // GIMP and CPDI require regular
    //  and qCPDI not allowed in 3D
	if(ElementBase::useGimp != POINT_GIMP)
    {   // using a GIMP method
		if(!mpmgrid.IsStructuredEqualElementsGrid())
		{	throw CommonException("GIMP needs a generated structured grid with equally sized elements",
								  "NairnMPM::ValidateOptions");
		}
        if(ElementBase::useGimp == QUADRATIC_CPDI)
        {   if(IsThreeD())
                throw CommonException("3D does not allow qCPDI shape functions; use lCPDI instead","NairnMPM::ValidateOptions");
        }
	}
    else
    {   // in Classic MPM or POINT_GIMP, cannot use traction BCs or +PS
        if(firstTractionPt!=NULL || firstFluxPt!=NULL || firstHeatFluxPt!=NULL)
			throw CommonException("Traction and flux boundary conditions require use of a GIMP MPM method.","NairnMPM::ValidateOptions");
		if(plusParticleSpin)
			throw CommonException("MPM+PS requires GIMP or CPDI","NairnMPM::ValidateOptions");
    }
    
    // Imperfect interface requires cartesian grid
	if(contact.hasImperfectInterface)
	{	if(!mpmgrid.IsStructuredEqualElementsGrid())
		{	throw CommonException("Imperfect interfaces need a generated structured grid with equally sized elements",
								  "NairnMPM::ValidateOptions");
		}
	}
	
    // 3D requires orthogonal grid and 1,8, or 27 particles per element
    // 2D requires 1,4, 9, 16, or 25 particles per element
	if(IsThreeD())
	{	if(mpmgrid.GetCartesian()!=CUBIC_GRID && mpmgrid.GetCartesian()!=ORTHOGONAL_GRID && mpmgrid.GetCartesian()!=VARIABLE_ORTHOGONAL_GRID)
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
            {   throw CommonException("Axisymmetric with Classic MPM needs an orthogonal grid","NairnMPM::ValidateOptions");
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
	{	if(!mpmgrid.IsStructuredEqualElementsGrid())
		{	// equally size elements needed for contact
			throw CommonException("Multimaterial mode needs a generated structured grid with equally sized elements","NairnMPM::ValidateOptions");
		}
	}
	
	// crack contact needs regular grid
	if(firstCrack!=NULL)
	{	if(!mpmgrid.IsStructuredEqualElementsGrid())
		{	// equally size elements needed for contact
			throw CommonException("Cracks need a generated structured grid with equally sized elements for contact calculations","NairnMPM::ValidateOptions");
		}
	}
	
	// The input commands have options to enter arbitray elements and nodes for an unstructured grid, but code
	// does nut current work with that style so it is isabled here. To use unstructure grid need at least:
	// 1. FindElementFromPoint() would need to search all elements
	// 2. Parallel patching would need rewriting
	if(!mpmgrid.IsStructuredGrid())
		throw CommonException("This code currently requies use of a generated structured grid","NairnMPM::ValidateOptions");
	
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
    return "NairnMPM";
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
            "    -H          Show this help\n"
            "    -np #       Set number of processors for parallel code\n"
            "    -r          Reverse byte order in archive files\n"
            "                   (default is to not reverse the bytes)\n"
            "    -v          Validate input file if DTD is provided in !DOCTYPE\n"
            "                   (default is to skip validation)\n"
            "    -w          Output results to current working directory\n"
            "                   (default is output to folder of input file)\n"
            "    -?          Show this help\n\n"
            "See http://osupdocs.forestry.oregonstate.edu/index.php/Main_Page for documentation\n\n"
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

// string about MPM additions (when printing out analysis type
const char *NairnMPM::MPMAugmentation(void) { return plusParticleSpin ? " +PS" : ""; }

