/*********************************************************************
    NairnMPM.cpp
    Nairn Research Group MPM Code
    
    Created by jnairn on Mon Nov 19 2001.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
*********************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "System/ArchiveData.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/CustomTask.hpp"
#include "Custom_Tasks/CalcJKTask.hpp"
#include "Custom_Tasks/PropagateTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "Exceptions/MPMTermination.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Boundary_Conditions/NodalTempBC.hpp"
#include "Boundary_Conditions/NodalConcBC.hpp"
#include "Boundary_Conditions/MatPtLoadBC.hpp"
#include "Elements/ElementBase.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Materials/RigidMaterial.hpp"
#include "Exceptions/MPMWarnings.hpp"
#include "Cracks/CrackNode.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include <time.h>

// global analysis object
NairnMPM *fmobj=NULL;

// to ignore crack interactions (only valid if 1 crack or non-interacting cracks)
//#define IGNORE_CRACK_INTERACTIONS

// global variables
double timestep=1.;			// time per MPM step (sec)
double strainTimestep;		// half time step when in US_AVG method
double mtime=0.;			// time for current step (sec)
double propTime=1e15;		// time interval between propagation calculations (sec)
int maxCrackFields=1;		// Maximum crack velocity fields at a node in a material (it is MAX_FIELDS_FOR_CRACKS if there are cracks)
int maxMaterialFields;		// Maximum velocity fields or number of independent materials

#pragma mark CONSTRUCTORS

// Contructor
NairnMPM::NairnMPM()
{
	version=7;						// main version
	subversion=2;					// subversion (must be < 10)
	buildnumber=1;					// build number
	mpmApproach=USAVG_METHOD;		// mpm method
	ptsPerElement=4;				// number of points per element (2D default, 3D changes it to 8)
	propagate=NO_PROPAGATION;				// default crack propagation type
	propagateDirection=DEFAULT_DIRECTION;	// default crack propagation direction
	propagateMat=0;							// default is new crack with no traction law
	hasTractionCracks=FALSE;		// if any crack segment has a traction law material
	maxtime=1.;						// maximum time (sec)
	FractCellTime=.5;				// fraction cell crossed in 1 step at wave spd
	mstep=0;						// step number
	volumeExtrap=FALSE;				// set if need volume extrapolations
	warnParticleLeftGrid=1;			// abort when this many leave the grid
	multiMaterialMode=false;		// multi-material mode
	
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

/*********************************************************************
    Main entry to read file and decode into objects
*********************************************************************/

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
	if(firstCrack!=NULL)
	{	if(propagate || archiver->WillArchiveJK(FALSE))
		{	nextTask=new CalcJKTask();
			if(propagate)
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
		
		sprintf(fline,"CPU Time per Step: %.3lf ms\n",1000.*cpuTime/((double)mstep));
		cout << fline;
	}
    
    //---------------------------------------------------
    // Trailer
    cout << "\n***** NairnMPM RUN COMPLETED\n";
}

/**********************************************************
	Main analysis loop for MPM analysis
**********************************************************/

void NairnMPM::MPMStep(void)
{
    long i,p,mi;
	int iel,numnds,matfld,nds[MaxShapeNds];
    short vfld;
    double mp,fn[MaxShapeNds],wt,xfrc,yfrc,xDeriv[MaxShapeNds],yDeriv[MaxShapeNds],zDeriv[MaxShapeNds];
    CrackHeader *nextCrack;
	TransportTask *nextTransport;
    char errMsg[100];
	MaterialBase *matID;
	
#pragma mark --- TASK 0: INITIALIZE
#ifdef LOG_PROGRESS
	char logLine[200];
	archiver->ClearLogFile();
	sprintf(logLine,"BEGIN STEP: %ld",mstep);
	archiver->WriteLogFile(logLine,NULL);
#endif
	
	// crack locations
	CrackField cfld[2];
	cfld[0].loc=NO_CRACK;		// NO_CRACK, ABOVE_CRACK, or BELOW_CRACK
	cfld[1].loc=NO_CRACK;
	
	// Zero Mass Matrix and vectors
	warnings.BeginStep();
	NodalPoint::ZeroAllNodesTask0();
	for(i=0;i<MaxShapeNds;i++) zDeriv[i]=0.;
    
    // Update forces applied to particles
	MatPtLoadBC::SetParticleFext(mtime);
	
	// undo dynamic velocity, temp, and conc BCs from rigid materials
	UnsetRigidBCs((BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
					(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
	UnsetRigidBCs((BoundaryCondition **)&firstTempBC,(BoundaryCondition **)&lastTempBC,
					(BoundaryCondition **)&firstRigidTempBC,(BoundaryCondition **)&reuseRigidTempBC);
	UnsetRigidBCs((BoundaryCondition **)&firstConcBC,(BoundaryCondition **)&lastConcBC,
					(BoundaryCondition **)&firstRigidConcBC,(BoundaryCondition **)&reuseRigidConcBC);
	
	// remove contact conditions
	CrackNode::RemoveCrackNodes();
        
    // turn off isothermal ramp when done and ramp step initialization
	thermal.CheckDone(mtime);
    
#pragma mark --- TASK 1: MASS MATRIX
#ifdef LOG_PROGRESS
	archiver->WriteLogFile("TASK 1: MASS MATRIX",NULL);
#endif

    /* Get mass matrix, find dimensionless particle locations,
            and find grid momenta
    */
    for(p=0;p<nmpms;p++)
	{	iel=mpm[p]->ElemID();
		
		// normal materials
		matID=theMaterials[mpm[p]->MatID()];
		if(!matID->Rigid())
		{	mp=mpm[p]->mp;			// in g
			matfld=matID->GetField();
		
			// get nodes and shape function for material point p
			if(multiMaterialMode)
				theElements[iel]->GetShapeFunctionsAndGradients(&numnds,fn,nds,&mpm[p]->pos,mpm[p]->GetNcpos(),xDeriv,yDeriv,zDeriv);
			else
				theElements[iel]->GetShapeFunctions(&numnds,fn,nds,&mpm[p]->pos,mpm[p]->GetNcpos());
			
			// get deformed particle volume if it will be needed (for transport tasks)
			if(volumeExtrap) mpm[p]->SetDilatedVolume();
			
			// Add particle property to each node in the element
			for(i=1;i<=numnds;i++)
			{	// Look for crack crossing and save until later
				if(firstCrack!=NULL)
				{	int cfound=0;
					Vector norm;
					cfld[0].loc=NO_CRACK;			// NO_CRACK, ABOVE_CRACK, or BELOW_CRACK
					cfld[1].loc=NO_CRACK;
					nextCrack=firstCrack;
					while(nextCrack!=NULL)
					{   vfld=nextCrack->CrackCross(mpm[p]->pos.x,mpm[p]->pos.y,nd[nds[i]]->x,nd[nds[i]]->y,&norm);
						if(vfld!=NO_CRACK)
						{	cfld[cfound].loc=vfld;
							cfld[cfound].norm=norm;
#ifdef IGNORE_CRACK_INTERACTIONS
							cfld[cfound].crackNum=1;	// appears to always be same crack, and stop when found one
							break;
#else
							cfld[cfound].crackNum=nextCrack->GetNumber();
							cfound++;
							if(cfound>1) break;			// stop if found two, if there are more then two, physics will be off
#endif
						}
						nextCrack=(CrackHeader *)nextCrack->GetNextObject();
					}
				}
				
				// momentum vector (and allocate velocity field if needed)
				vfld=nd[nds[i]]->AddMomentumTask1(matfld,cfld,fn[i]*mp,&mpm[p]->vel);
				mpm[p]->vfld[i]=vfld;
				
				// add to lumped mass matrix
				nd[nds[i]]->AddMassTask1(vfld,matfld,mp*fn[i]);
				
				// crack contact calculations
				contact.AddDisplacementVolumeTask1(vfld,matfld,nd[nds[i]],mpm[p],fn[i]);
				
				// transport calculations
				nextTransport=transportTasks;
				while(nextTransport!=NULL)
					nextTransport=nextTransport->Task1Extrapolation(nd[nds[i]],mpm[p],fn[i]);
				
				// material contact calculations
				if(multiMaterialMode)
					nd[nds[i]]->AddMassGradient(vfld,matfld,mp,xDeriv[i],yDeriv[i],zDeriv[i],mpm[p]);
			}
		}
		
		// For Rigid materials create velocity BC on each node in the element
		else
		{	numnds=theElements[iel]->NumberNodes();
			double rvalue;
			for(i=1;i<=numnds;i++)
			{   mi=theElements[iel]->nodes[i-1];		// 1 based node
				RigidMaterial *rigid=(RigidMaterial *)(theMaterials[mpm[p]->MatID()]);
				
				// check skewed, x or y direction velocities
				if(rigid->RigidDirection(X_DIRECTION+Y_DIRECTION))
				{	// get magnitude and angle (in radians, cw from x+ axis)
					double vel=sqrt(mpm[p]->vel.x*mpm[p]->vel.x+mpm[p]->vel.y*mpm[p]->vel.y);
					double angle;
					if(DbleEqual(mpm[p]->vel.x,0.))
						angle = (mpm[p]->vel.y>0) ? -PI_CONSTANT/2. : PI_CONSTANT/2. ;
					else if(mpm[p]->vel.x>0)
						angle=-atan(mpm[p]->vel.y/mpm[p]->vel.x);
					else
						angle=PI_CONSTANT-atan(mpm[p]->vel.y/mpm[p]->vel.x);
					
					// adjust magnitude if has setting function and not currently zero
					// but setting function for skewed conditions must be positive otherwise skew direction
					// toggles back and forth, rather than continuing in same directioon
					if(rigid->GetSetting(&rvalue,mtime))
					{	if(!DbleEqual(vel,0.))
						{	vel=rvalue;
							mpm[p]->vel.x=vel*cos(angle);
							mpm[p]->vel.y=-vel*sin(angle);
						}
					}
					
					if(DbleEqual(mpm[p]->vel.x,0.) || DbleEqual(mpm[p]->vel.y,0.))
					{	// special case of on-axis zero fixes both directions
						SetRigidBCs(mi,X_DIRECTION,mpm[p]->vel.x,0.,
								(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
								(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
						SetRigidBCs(mi,Y_DIRECTION,mpm[p]->vel.y,0.,
								(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
								(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
					}
					else
					{	angle*=180./PI_CONSTANT;			// convert to degrees for BC methods
						SetRigidBCs(mi,SKEW_DIRECTION,vel,angle,
								(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
								(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
					}
				}
				else if(rigid->RigidDirection(X_DIRECTION))
				{	if(rigid->GetSetting(&rvalue,mtime)) mpm[p]->vel.x=rvalue;
					SetRigidBCs(mi,X_DIRECTION,mpm[p]->vel.x,0.,
							(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
							(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
				}
				else if(rigid->RigidDirection(Y_DIRECTION))
				{	if(rigid->GetSetting(&rvalue,mtime)) mpm[p]->vel.y=rvalue;
					SetRigidBCs(mi,Y_DIRECTION,mpm[p]->vel.y,0.,
							(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
							(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
				}
				
				// z direction
				if(rigid->RigidDirection(Z_DIRECTION))
				{	if(rigid->GetSetting(&rvalue,mtime)) mpm[p]->vel.z=rvalue;
					SetRigidBCs(mi,Z_DIRECTION,mpm[p]->vel.z,0.,
							(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
							(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
				}
					
				// temperature
				if(rigid->RigidTemperature())
				{	if(rigid->GetSetting(&rvalue,mtime)) mpm[p]->pTemperature=rvalue;
					SetRigidBCs(mi,TEMP_DIRECTION,mpm[p]->pTemperature,0.,
							(BoundaryCondition **)&firstTempBC,(BoundaryCondition **)&lastTempBC,
							(BoundaryCondition **)&firstRigidTempBC,(BoundaryCondition **)&reuseRigidTempBC);
				}
				
				// concentration
				if(rigid->RigidConcentration())
				{	if(rigid->GetSetting(&rvalue,mtime)) mpm[p]->pConcentration=rvalue;
					SetRigidBCs(mi,CONC_DIRECTION,mpm[p]->pConcentration,0.,
							(BoundaryCondition **)&firstConcBC,(BoundaryCondition **)&lastConcBC,
							(BoundaryCondition **)&firstRigidConcBC,(BoundaryCondition **)&reuseRigidConcBC);
				}
			}
		}
    }
	
	// if any left over rigid BCs, delete them now
	RemoveRigidBCs((BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,(BoundaryCondition **)&firstRigidVelocityBC);
	RemoveRigidBCs((BoundaryCondition **)&firstTempBC,(BoundaryCondition **)&lastTempBC,(BoundaryCondition **)&firstRigidTempBC);
	RemoveRigidBCs((BoundaryCondition **)&firstConcBC,(BoundaryCondition **)&lastConcBC,(BoundaryCondition **)&firstRigidConcBC);
	
	// total nodal masses and clear out unused fields - in case needed, e.g. for cracks
	NodalPoint::GetNodalMassesTask1();
	
#pragma mark --- TASK 1b: TRANSPORT EXTRAS
#ifdef LOG_PROGRESS
	archiver->WriteLogFile("TASK 1b: TRANSPORT EXTRAS",NULL);
#endif

	// Extra calculations for transport
	nextTransport=transportTasks;
	while(nextTransport!=NULL)
		nextTransport=nextTransport->GetValuesAndGradients(mtime);
	
#pragma mark --- TASK 2: UPDATE STRAINS
#ifdef LOG_PROGRESS
	archiver->WriteLogFile("TASK 2: UPDATE STRAINS",NULL);
#endif

	/* Adjust extrapolated velocity for contact, impose velocity BCs, and then
		update strains with those velocities
		NOTE: Switch order of contact and BCs (8/12/2009)
    */
	NodalPoint::MaterialContact(multiMaterialMode,FALSE,timestep);
	CrackHeader::ContactConditions(TRUE);
	NodalVelBC::GridMomentumConditions(TRUE);
    if(mpmApproach==USF_METHOD || mpmApproach==USAVG_METHOD)
        MPMBase::FullStrainUpdate(strainTimestep,FALSE,np);
	
#pragma mark --- TASK 3: FORCES
#ifdef LOG_PROGRESS
	archiver->WriteLogFile("TASK 3: FORCES",NULL);
#endif
	
	// Get total grid point forces (except external forces)
    for(p=0;p<nmpms;p++)
	{	matID=theMaterials[mpm[p]->MatID()];
		if(matID->Rigid()) continue;
	
		// get transport tensors (if needed)
		if(transportTasks!=NULL)
			matID->LoadTransportProps(mpm[p]);
			
        mp=mpm[p]->mp;					// in g
		matfld=matID->GetField();		// material field
		
        // find shape functions and derviatives
    	iel=mpm[p]->ElemID();
		theElements[iel]->GetShapeGradients(&numnds,fn,nds,mpm[p]->GetNcpos(),xDeriv,yDeriv,zDeriv);

        // Add particle property to each node in the element
        for(i=1;i<=numnds;i++)
		{	vfld=(short)mpm[p]->vfld[i];					// crack velocity field to use
    
            /* internal force vector (in g mm/sec^2)
                            (note: stress is specific stress in units N/m^2 cm^3/g
                                            Multiply by 1000 to make it mm/sec^2)
            */
			Vector theFint=mpm[p]->Fint(xDeriv[i],yDeriv[i],zDeriv[i]);
			nd[nds[i]]->AddFintTask3(vfld,matfld,mpm[p]->Fint(xDeriv[i],yDeriv[i],zDeriv[i]));
            
            // body forces (not 3D yet)
			if(bodyFrc.GetGravity(&xfrc,&yfrc))
				nd[nds[i]]->AddFintTask3(vfld,matfld,MakeVector(mp*fn[i]*xfrc,mp*fn[i]*yfrc,0.));
            
			// external force vector
            nd[nds[i]]->AddFextTask3(vfld,matfld,mpm[p]->Fext(fn[i]));
			
			// transport forces
			nextTransport=transportTasks;
			while(nextTransport!=NULL)
				nextTransport=nextTransport->AddForces(nd[nds[i]],mpm[p],fn[i],xDeriv[i],yDeriv[i],zDeriv[i]);
        }
		
		// clear coupled dissipated energy
		if(ConductionTask::energyCoupling) mpm[p]->SetDispEnergy(0.);
    }
	
	// traction law forces
	if(hasTractionCracks)
	{	nextCrack=firstCrack;
		while(nextCrack!=NULL)
		{	nextCrack->TractionFext();
			nextCrack=(CrackHeader *)nextCrack->GetNextObject();
		}
	}
	
	// crack tip heating
	if(conduction) conduction->AddCrackTipHeating();
	
	// interface forces
	CrackNode::InterfaceOnKnownNodes();
    
	// Find Grid total force
	NodalPoint::GetGridForcesTask3(bodyFrc.GetDamping());
	
    /* Imposed BCs on ftot to get correct grid BCs for velocity
		and concentration and temperature.
	*/
    NodalVelBC::ConsistentGridForces();
	nextTransport=transportTasks;
	while(nextTransport!=NULL)
		nextTransport=nextTransport->SetTransportForceBCs(timestep);
    
#pragma mark --- TASK 4: UPDATE MOMENTA
#ifdef LOG_PROGRESS
	archiver->WriteLogFile("TASK 4: UPDATE MOMENTA",NULL);
#endif

	// Update grid momenta and transport rates
	NodalPoint::UpdateGridMomentaTask4(timestep);
	NodalPoint::MaterialContact(multiMaterialMode,TRUE,timestep);
	CrackNode::CrackContactTask4(timestep);
	nextTransport=transportTasks;
	while(nextTransport!=NULL)
		nextTransport=nextTransport->TransportRates(timestep);
	
#pragma mark --- TASK 5: UPDATE PARTICLES
#ifdef LOG_PROGRESS
	archiver->WriteLogFile("TASK 5: UPDATE PARTICLES",NULL);
#endif


    // Update particle position, velocity, temp, and conc
	Vector delv,*acc;
	bodyFrc.TrackAlpha();
    for(p=0;p<nmpms;p++)
	{	matID=theMaterials[mpm[p]->MatID()];
		if(!matID->Rigid())
		{	// get shape functions
			iel=mpm[p]->ElemID();
			theElements[iel]->GetShapeFunctions(&numnds,fn,nds,mpm[p]->GetNcpos());

			// Update particle position and velocity
			matfld=matID->GetField();
			acc=mpm[p]->GetAcc();
			ZeroVector(acc);
			ZeroVector(&delv);
			nextTransport=transportTasks;
			while(nextTransport!=NULL)
				nextTransport=nextTransport->ZeroTransportRate();
			for(i=1;i<=numnds;i++)
			{	nd[nds[i]]->IncrementDelvaTask5((short)mpm[p]->vfld[i],matfld,fn[i],&delv,acc);
				nextTransport=transportTasks;
				while(nextTransport!=NULL)
					nextTransport=nextTransport->IncrementTransportRate(nd[nds[i]],fn[i]);
			}
						
			// update position in mm and velocity in mm/sec
			mpm[p]->MovePosition(timestep,&delv);
			mpm[p]->MoveVelocity(timestep,bodyFrc.GetAlpha(),&delv);
			
			// update transport values
			nextTransport=transportTasks;
			while(nextTransport!=NULL)
				nextTransport=nextTransport->MoveTransportValue(mpm[p],timestep);
			
			// thermal ramp
			thermal.UpdateParticleTemperature(&mpm[p]->pTemperature,timestep);
			
			// update feedback coefficient
			bodyFrc.TrackAlpha(&mpm[p]->vel,mpm[p]->mp);
		}
		
		else
		{	// rigid materials at constant velocity
			mpm[p]->MovePosition(timestep,&mpm[p]->vel);
		}
	}
	
	// update damping coefficient
	bodyFrc.UpdateAlpha(timestep,mtime);
	
#pragma mark --- TASK 6: UPDATE STRAINS
#ifdef LOG_PROGRESS
	archiver->WriteLogFile("TASK 6: UPDATE STRAINS",NULL);
#endif

    /* For SZS and USAVG Methods
            Extrapolate new particle velocities to grid
            Update stress and strain
    */
    if(mpmApproach==SZS_METHOD || mpmApproach==USAVG_METHOD)
    {	// zero again
		NodalPoint::RezeroAllNodesTask6();
        
        // loop over particles
        for(p=0;p<nmpms;p++)
		{	matID=theMaterials[mpm[p]->MatID()];
			if(matID->Rigid()) continue;
            mp=mpm[p]->mp;			// in g
			matfld=matID->GetField();
    
            // find shape functions
			iel=mpm[p]->ElemID();
			if(multiMaterialMode)
				theElements[iel]->GetShapeGradients(&numnds,fn,nds,mpm[p]->GetNcpos(),xDeriv,yDeriv,zDeriv);
			else
				theElements[iel]->GetShapeFunctions(&numnds,fn,nds,mpm[p]->GetNcpos());
            
            // update vector for velocity extrapolation
            for(i=1;i<=numnds;i++)
			{	vfld=(short)mpm[p]->vfld[i];
				
            	// velocity from updated velocities
                nd[nds[i]]->AddMomentumTask6(vfld,matfld,fn[i]*mp,&mpm[p]->vel);
				
				// add updated displacement (if cracks, not 3D)
				contact.AddDisplacementTask6(vfld,matfld,nd[nds[i]],mpm[p],fn[i]);
				
				// material contact calculations
				if(multiMaterialMode)
					nd[nds[i]]->AddMassGradient(vfld,matfld,mp,xDeriv[i],yDeriv[i],zDeriv[i],mpm[p]);
            }
        }
		
		// update nodal values for transport properties (when coupled to strain)
		nextTransport=transportTasks;
		while(nextTransport!=NULL)
			nextTransport=nextTransport->UpdateNodalValues(timestep);
        
        // Update strains with newly extrapolated momenta
		NodalPoint::MaterialContact(multiMaterialMode,FALSE,timestep);
		CrackHeader::ContactConditions(FALSE);
		NodalVelBC::GridMomentumConditions(FALSE);
        MPMBase::FullStrainUpdate(strainTimestep,(mpmApproach==USAVG_METHOD),np);
    }
   
#pragma mark --- TASK 7a: CUSTOM TASKS START
#ifdef LOG_PROGRESS
	archiver->WriteLogFile("TASK 7a: CUSTOM TASKS START",NULL);
#endif

    /* Call all tasks. The tasks can do initializes needed
            for this step. If any task needs nodal extrapolations
            they should set needNodalExtraps to TRUE. If it does
			not need them, leave it alone
    */
    bool needExtrapolations=FALSE;
    CustomTask *nextTask=theTasks;
    while(nextTask!=NULL)
	{	bool taskNeedsExtrapolations=FALSE;
    	nextTask=nextTask->PrepareForStep(taskNeedsExtrapolations);
		// if it was set to TRUE, trasfer to to global setting
		if(taskNeedsExtrapolations) needExtrapolations=TRUE;
	}
        
#pragma mark --- TASK 7b: CUSTOM EXTRAPOLATIONS
#ifdef LOG_PROGRESS
	archiver->WriteLogFile("TASK 7b: CUSTOM EXTRAPOLATIONS",NULL);
#endif

    /* Extrapolate particle info to grid if needed for
                    any custom task
    */
    if(needExtrapolations)
    {	// call each task for initialization prior to extrapolations
    	nextTask=theTasks;
        while(nextTask!=NULL)
            nextTask=nextTask->BeginExtrapolations();
        
        // particle loop
        for(p=0;p<nmpms;p++)
        {   // Load element coordinates
			matID=theMaterials[mpm[p]->MatID()];
			if(matID->Rigid()) continue;
			matfld=matID->GetField();
    
            // find shape functions and derviatives
			iel=mpm[p]->ElemID();
			theElements[iel]->GetShapeGradients(&numnds,fn,nds,mpm[p]->GetNcpos(),xDeriv,yDeriv,zDeriv);
            
			// Add particle property to each node in the element
            for(i=1;i<=numnds;i++)
            {   // global mass matrix
				vfld=(short)mpm[p]->vfld[i];				// velocity field to use
                wt=fn[i]*mpm[p]->mp;
                
                // possible extrapolation to the nodes
                nextTask=theTasks;
                while(nextTask!=NULL)
                    nextTask=nextTask->NodalExtrapolation(nd[nds[i]],mpm[p],vfld,matfld,wt);
                    
                // possible extrapolation to the particle
                nextTask=theTasks;
                while(nextTask!=NULL)
                    nextTask=nextTask->ParticleCalculation(nd[nds[i]],mpm[p],vfld,matfld,fn[i],xDeriv[i],yDeriv[i],zDeriv[i]);
            }
            
            // possible single calculations for each particle
            nextTask=theTasks;
            while(nextTask!=NULL)
                nextTask=nextTask->ParticleExtrapolation(mpm[p]);
        }
        
        // finished with extrapolations
        nextTask=theTasks;
        while(nextTask!=NULL)
            nextTask=nextTask->EndExtrapolations();
    }
        
#pragma mark --- TASK 7c: CUSTOM CALCULATIONS
#ifdef LOG_PROGRESS
	archiver->WriteLogFile("TASK 7c: CUSTOM CALCULATIONS",NULL);
#endif

    // Do the custom calculations
    nextTask=theTasks;
    while(nextTask!=NULL)
        nextTask=nextTask->StepCalculation();
        
#pragma mark --- TASK 7d: CUSTOM FINISH
#ifdef LOG_PROGRESS
	archiver->WriteLogFile("TASK 7d: CUSTOM FINISH",NULL);
#endif

    // Call tasks in case any need to clean up
    nextTask=theTasks;
    while(nextTask!=NULL)
    	nextTask=nextTask->FinishForStep();
        
#pragma mark --- TASK 8a: MOVE CRACKS
#ifdef LOG_PROGRESS
	archiver->WriteLogFile("TASK 8a: MOVE CRACKS",NULL);
#endif

	// Move crack surfaces
    if(firstCrack!=NULL)
	{	// prepare multimaterial fields for moving cracks
		
    	nextCrack=firstCrack;
    	while(nextCrack!=NULL)
        {   if(!nextCrack->MoveCrack(ABOVE_CRACK))
			{	sprintf(errMsg,"Crack No. %d surface (above) has left the grid.",nextCrack->GetNumber());
            	throw MPMTermination(errMsg,"NairnMPM::MPMStep");
			}
            if(!nextCrack->MoveCrack(BELOW_CRACK))
			{	sprintf(errMsg,"Crack No. %d surface (below) has left the grid.",nextCrack->GetNumber());
            	throw MPMTermination(errMsg,"NairnMPM::MPMStep");
			}
            nextCrack=(CrackHeader *)nextCrack->GetNextObject();
        }
		
		if(!contact.GetMoveOnlySurfaces()) NodalPoint::GetGridCMVelocitiesTask8();
        nextCrack=firstCrack;
        while(nextCrack!=NULL)
        {   if(!nextCrack->MoveCrack())
			{	sprintf(errMsg,"Crack No. %d position or surface has left the grid.",nextCrack->GetNumber());
            	throw MPMTermination(errMsg,"NairnMPM::MPMStep");
			}
            nextCrack=(CrackHeader *)nextCrack->GetNextObject();
        }
		
		// update crack tractions
		if(hasTractionCracks)
		{	nextCrack=firstCrack;
			while(nextCrack!=NULL)
			{	nextCrack->UpdateCrackTractions();
				nextCrack=(CrackHeader *)nextCrack->GetNextObject();
			}
		}
    }
	
#pragma mark --- TASK 8b: RESET ELEMEMTS
#ifdef LOG_PROGRESS
	archiver->WriteLogFile("TASK 8b: RESET ELEMEMTS",NULL);
#endif

    // See if any particles have changed elements
	// Stop if off the grid
    for(p=0;p<nmpms;p++)
    {	if(!ResetElement(mpm[p]))
		{	if(warnings.Issue(warnParticleLeftGrid,-1)==REACHED_MAX_WARNINGS)
			{   sprintf(errMsg,"Particle No. %ld left the grid\n  (plot x displacement to see it).",p+1);
				mpm[p]->origpos.x=-1.e6;
				throw MPMTermination(errMsg,"NairnMPM::MPMStep");
			}
			
			// bring back to the previous element
			ReturnToElement(mpm[p]);
        }
    }
}

/**********************************************************
    Set boundary conditions determined by moving
	rigid paticles
**********************************************************/

void NairnMPM::SetRigidBCs(long mi,int type,double value,double angle,BoundaryCondition **firstBC,
						BoundaryCondition **lastBC,BoundaryCondition **firstRigidBC,BoundaryCondition **reuseRigidBC)
{
	BoundaryCondition *newBC=NULL;
	NodalVelBC *velBC;
	
	// check if already set in that direction by actual BC or by previous rigid BC
	// New rigid BC's can only be on free directions
	if(nd[mi]->fixedDirection&type) return;
	
	// create new boundary conditions
	switch(type)
	{	case SKEW_DIRECTION:
			if(nd[mi]->fixedDirection&(X_DIRECTION+Y_DIRECTION)) return;
			if(*reuseRigidBC!=NULL)
				velBC=(NodalVelBC *)((*reuseRigidBC)->SetRigidProperties(mi,type,CONSTANT_VALUE,value));
			else
			{	velBC=new NodalVelBC(mi,type,CONSTANT_VALUE,value,(double)0.);
				if(velBC==NULL) throw CommonException("Memory error allocating rigid particle boundary condition.",
													  "NairnMPM::SetRigidBCs");
			}
			velBC->SetSkewAngle(angle);
			newBC=(BoundaryCondition *)velBC;
			break;
	
		case X_DIRECTION:
		case Y_DIRECTION:
		case Z_DIRECTION:
			if(*reuseRigidBC!=NULL)
				newBC=(*reuseRigidBC)->SetRigidProperties(mi,type,CONSTANT_VALUE,value);
			else
			{	newBC=(BoundaryCondition *)(new NodalVelBC(mi,type,CONSTANT_VALUE,value,(double)0.));
				if(newBC==NULL) throw CommonException("Memory error allocating rigid particle boundary condition.",
													  "NairnMPM::SetRigidBCs");
			}
			break;
			
		case TEMP_DIRECTION:
			if(*reuseRigidBC!=NULL)
				newBC=(*reuseRigidBC)->SetRigidProperties(mi,type,CONSTANT_VALUE,value);
			else
			{	newBC=(BoundaryCondition *)(new NodalTempBC(mi,CONSTANT_VALUE,value,(double)0.));
				if(newBC==NULL) throw CommonException("Memory error allocating rigid particle boundary condition.",
													  "NairnMPM::SetRigidBCs");
			}
			break;
			
		case CONC_DIRECTION:
			if(*reuseRigidBC!=NULL)
				newBC=(*reuseRigidBC)->SetRigidProperties(mi,type,CONSTANT_VALUE,value);
			else
			{	newBC=(BoundaryCondition *)(new NodalConcBC(mi,CONSTANT_VALUE,value,(double)0.));
				if(newBC==NULL) throw CommonException("Memory error allocating rigid particle boundary condition.",
													  "NairnMPM::SetRigidBCs");
			}
			break;
			
		default:
			break;
	}
	
	// *firstBC and *lastBC are first and last of this type
	// *firstRigidBC will save the first one
	// if *reuseRigidBC!=NULL, then reusing previous rigid BCs
	if(*firstBC==NULL)
	{	// Only happens when no normal BCs and no rigidBCs to reuse so start with this one
		*firstBC=newBC;
		*firstRigidBC=newBC;
		// reuseRigidBC must be NULL already
	}
	else
	{	if(*reuseRigidBC!=NULL)
		{	// next object of last BC is already set
			// firstRigidBC is already valid
			// advance to reuse next one (or could get to NULL if all used up)
			*reuseRigidBC=(BoundaryCondition *)(*reuseRigidBC)->GetNextObject();
		}
		else
		{	// created a new one or ran out of ones to reuse
			(*lastBC)->SetNextObject(newBC);
			if(*firstRigidBC==NULL) *firstRigidBC=newBC;
		}
	}
	*lastBC=newBC;
}

/**********************************************************
	Unset nodal dof for rigid BCs so can try to reuse
	them without needing new memory allocations
**********************************************************/

void NairnMPM::UnsetRigidBCs(BoundaryCondition **firstBC,BoundaryCondition **lastBC,
									BoundaryCondition **firstRigidBC,BoundaryCondition **reuseRigidBC)
{
	// exit if none
	if(*firstRigidBC==NULL) return;
	
	// unset dynamic ones
    BoundaryCondition *nextBC=*firstRigidBC;
	while(nextBC!=NULL)
		nextBC=nextBC->UnsetDirection();
	
	// were they all rigid BCs, but still has some?
	if(*firstBC==*firstRigidBC)
	{	*lastBC=NULL;				// since none set yet
		*reuseRigidBC=*firstRigidBC;
		return;
	}
	
	// otherwise search for last actual grid BC
    nextBC=*firstBC;
	while(nextBC->GetNextObject()!=*firstRigidBC)
		nextBC=(BoundaryCondition *)nextBC->GetNextObject();
	*lastBC=nextBC;
	*reuseRigidBC=*firstRigidBC;
}

/**********************************************************
	Debugging aid to count BCs on each step
**********************************************************/

void NairnMPM::CountBCs(BoundaryCondition **firstBC,BoundaryCondition **lastBC,BoundaryCondition **firstRigidBC)
{
	if(*firstBC==NULL) return;
	
	BoundaryCondition *nextBC=*firstBC;
	int count=0,fullcount=0;
	while(nextBC!=NULL)
	{	if(nextBC==*firstRigidBC)
		{	cout << "# " << count << " on grid, ";
			count=0;
		}
		count++;
		fullcount++;
		if(nextBC==*lastBC)
		{	cout << count << " rigid, ";
			count=0;
		}
		nextBC=(BoundaryCondition *)nextBC->GetNextObject();
	}
	cout << count << " undeleted, total = " << fullcount << endl;
}
	
/**********************************************************
	Remove any dynamically created boundary conditions
	that are no longer needed
**********************************************************/

void NairnMPM::RemoveRigidBCs(BoundaryCondition **firstBC,BoundaryCondition **lastBC,BoundaryCondition **firstRigidBC)
{
	// exit if none
	if(*firstRigidBC==NULL) return;
	
	// trap if this step did not reuse any BCs
	BoundaryCondition *nextBC,*prevBC;
	if(*lastBC==NULL)
	{	nextBC=*firstBC;		// will be deleting them all below
		// None reused and no grid ones either?
		*firstBC=NULL;
		*firstRigidBC=NULL;
	}
	else
	{	nextBC=(BoundaryCondition *)(*lastBC)->GetNextObject();
		(*lastBC)->SetNextObject(NULL);			// next one on lastBC needs to be NULL
		// Were none of the rigid ones resued? (they are deleted below)
		if(nextBC==*firstRigidBC) *firstRigidBC=NULL;
	}
	
	// delete any rigid BCs that were not reused
	while(nextBC!=NULL)
	{	prevBC=nextBC;
		nextBC=(BoundaryCondition *)prevBC->GetNextObject();
		delete prevBC;
	}
}

/**********************************************************
	Find element for particle. Return FALSE if left
	the grid or for GIMP moved to an edge element
**********************************************************/

int NairnMPM::ResetElement(MPMBase *mpt)
{
    int i;
	
    // check current element
    if(theElements[mpt->ElemID()]->PtInElement(mpt->pos)) return TRUE;
    
    // check others
    for(i=0;i<nelems;i++)
    {	if(theElements[i]->PtInElement(mpt->pos))
		{	if(theElements[i]->OnTheEdge()) return FALSE;
			mpt->ChangeElemID(i);
            return TRUE;
        }
    }
    return FALSE;
}

/**********************************************************
	Find element for particle. Return FALSE if left
	the grid or for GIMP moved to an edge element
**********************************************************/

void NairnMPM::ReturnToElement(MPMBase *mpt)
{
	Vector outside=mpt->pos;
    int elemID=mpt->ElemID();
	Vector inside,middle,origin;
	int pass;
	
	// try to retrace position, if fails, take element centroid
	inside.x=outside.x-timestep*mpt->vel.x;
	inside.y=outside.y-timestep*mpt->vel.y;
	inside.z=outside.z-timestep*mpt->vel.z;
	if(!theElements[elemID]->PtInElement(inside))
		theElements[elemID]->GetXYZCentroid(&inside);
	origin=inside;
		
	// bisect 10 times
	for(pass=1;pass<=10;pass++)
	{	middle.x=(outside.x+inside.x)/2.;
		middle.y=(outside.y+inside.y)/2.;
		middle.z=(outside.z+inside.z)/2.;
		
		if(theElements[elemID]->PtInElement(middle))
			inside=middle;
		else
			outside=middle;
	}
	
	// move to inside
	mpt->SetPosition(&inside);
	
	// change velocity for movement form starting position to new edge position (but seems to not be good idea)
	//mpt->vel.x=(inside.x-origin.x)/timestep;
	//mpt->vel.y=(inside.y-origin.y)/timestep;
	//mpt->vel.z=(inside.z-origin.z)/timestep;
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
    long p,i;
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
		if(theMaterials[matid]->Rigid()) continue;
		maxMaterialFields=theMaterials[matid]->SetField(maxMaterialFields,multiMaterialMode,matid);
	
		// element and mp properties
		if(IsThreeD())
		{	volume=theElements[mpm[p]->ElemID()]->GetVolume()/1000.;	// in cm^2
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
        
        // check time step
        crot=theMaterials[matid]->WaveSpeed()/10.;		// in cm/sec
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
	if(maxMaterialFields==0)
		throw CommonException("No material points found with a non-rigid material","NairnMPM::PreliminaryCalcs");
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
    
    sprintf(fline,"Number of Material Points: %ld",nmpms);
    cout << fline << endl;
	
	// background grid info
	mpmgrid.Output(ptsPerElement);
    
    sprintf(fline,"Adjusted time step (ms): %.7e",1000.*timestep);
    cout << fline << endl;
    
	// contact law materials and cracks
	contact.SetNormalCODCutoff(mpmgrid.GetMinCellDimension());
	
    // progation time step and other settings when has cracks
    if(firstCrack!=NULL)
	{	if(propagate)
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
	warnParticleLeftGrid=warnings.CreateWarning("particle left the grid",(long)warnParticleLeftGrid);
	
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
    cout << "    Expects Xerces 2.8.0 or newer" << endl;
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

