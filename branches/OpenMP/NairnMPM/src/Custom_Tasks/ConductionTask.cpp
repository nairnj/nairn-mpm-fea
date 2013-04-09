/*********************************************************************************************
    ConductionTask.cpp
    NairnMPM
    
    Created by John Nairn on Fri Oct 15 2004
    Copyright (c) 2004 John A. Nairn, All rights reserved.
 
	Conduction calculations
   -------------------------
	Initialization:
		Set gTemperature, gMpCp, fcond on node to zero;
	Mass and Momentum Task
		Extrapolate gTemperature and gMpCp (Task1Extrapolation())
		Divide gTemperature by gMpCp (GetNodalValue()
		Impose grid T BCs (ImposeValueBCs())
		Find grad T on particles (GetGradients())
	Grid Forces Task
		Extrapolate conductivity force to fcond (AddForces())
			(include heat sources, and energy coupling)
		Add crack tip heating to fcond (AddCrackTipHeating())
		Finish grid BCs and impose flux BCs in fcond (SetTransportForceBCs())
	Update Momenta Task
		Divide fcond by gMpCp to get temperature rates (TransportRates())
	Update Particles Task
		Each particle: zero rate, extrapolate from nodes to particle, then
			update particle (IncrementTransportRate(),MoveTransportRate()).
	Update strains last (if used)
		Update gTemperature on nodes (UpdateNodalValues())
	Update strains on particles (both places)
		Extrapolate grid temperature to particle using IncrementValueExtrap()
		Find dT = extrapolated value minus previous extrapolated value and then
			store extrapolated value in pPreviousTemperature (GetDeltaValue())
		This task is done because contititutive laws work better when temperature comes
			from grid extrapolation instead of particle temperature. The current grid
			based result is stored in a different particle temperature. The particle
			temperature, however, is the one used in conduction calculations. Using the
			grid based one causes numerical diffision.
*********************************************************************************************/

#include "Custom_Tasks/ConductionTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "Boundary_Conditions/NodalTempBC.hpp"
#include "Cracks/CrackHeader.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Cracks/CrackSegment.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Materials/RigidMaterial.hpp"

// global
bool ConductionTask::active=FALSE;
bool ConductionTask::crackTipHeating=FALSE;
bool ConductionTask::energyCoupling=FALSE;
bool ConductionTask::AVHeating=TRUE;
ConductionTask *conduction=NULL;

#pragma mark INITIALIZE

// Constructors
ConductionTask::ConductionTask()
{	// allocate diffusion data on each particle
    int p;
	for(p=0;p<nmpmsNR;p++)
		mpm[p]->AllocateTemperature();
}

// Return name of this task
const char *ConductionTask::TaskName(void) { return "conduction analysis"; }

#pragma mark STANDARD METHODS

// conduction analysis settings
TransportTask *ConductionTask::TransportOutput(void)
{	TransportTask::TransportOutput();
	if(crackTipHeating)
		cout << "   Crack tip heating activated" << endl;
	return nextTask;
}

// conduction analysis settings
void ConductionTask::ThermodynamicsOutput(void)
{   // the system
    if(ConductionTask::IsSystemIsolated())
        cout << "System: isolated" << endl;
    else
        cout << "System: nonisolated" << endl;;
    if(active)
        cout << "Particles: nonisolated";
    else
        cout << "Particles: isolated";
    if(energyCoupling)
        cout << " and adiabatic" << endl;
    else
        cout << " and isothermal" << endl;
}

// is the system isolated?
// update this if new thermal BCs are added
bool ConductionTask::IsSystemIsolated(void)
{
    // if has ramp, then is it not isolated
    if(thermal.Active()) return FALSE;
    
    // if no conduction then is isolated
    if(!active) return TRUE;
    
    // if conduction is active, still isolated if no thermal BC
    if(firstTempBC!=NULL) return FALSE;
    
    // check if any active rigid particles set temperature
	if(RigidMaterial::someSetTemperature) return FALSE;
    
    // must be isolated
    return TRUE;
}


// adjust time for given cell size if needed
TransportTask *ConductionTask::TransportTimeStep(int matid,double dcell,double *tmin)
{	double diffCon=theMaterials[matid]->MaximumDiffusivity();
	double tst=(dcell*dcell)/(4.*diffCon);						// factor 2 shorter than minimum
	if(tst<*tmin) *tmin=tst;
	return nextTask;
}

#pragma mark MASS AND MOMENTUM EXTRAPOLATIONS

// Task 1 Extrapolation of concentration to the grid
TransportTask *ConductionTask::Task1Extrapolation(NodalPoint *ndpt,MPMBase *mptr,double shape)
{	double Cp=theMaterials[mptr->MatID()]->GetHeatCapacity(mptr);   // mJ/(g-K)
	double arg = mptr->mp*Cp*shape;                                 // mJ/K
	ndpt->gTemperature+=mptr->pTemperature*arg;                     // mJ
	ndpt->gMpCp+=arg;                                               // mJ/K
	return nextTask;
}

// Task 1b - get grid temperature on one node
TransportTask *ConductionTask::GetNodalValue(NodalPoint *ndptr)
{   if(ndptr->NumberNonrigidParticles()>0)
		ndptr->gTemperature /= ndptr->gMpCp;
	return nextTask;
}


// Task 1b - impose grid-based temperature BCs
void ConductionTask::ImposeValueBCs(double stepTime)
{
	int i;
		
	// Copy no-BC temperature
    NodalTempBC *nextBC=firstTempBC;
    while(nextBC!=NULL)
		nextBC=nextBC->CopyNodalTemperature(nd[nextBC->GetNodeNum()]);
	
    // Zero them all
	double mstime=1000.*stepTime;
    nextBC=firstTempBC;
    while(nextBC!=NULL)
	{   i=nextBC->GetNodeNum(mstime);
		if(i!=0) nd[i]->gTemperature=0.;
        nextBC=(NodalTempBC *)nextBC->GetNextObject();
    }
	
    // Now add all temperature to nodes with temperature BCs
    nextBC=firstTempBC;
    while(nextBC!=NULL)
	{   i=nextBC->GetNodeNum(mstime);
		if(i!=0) nd[i]->gTemperature+=nextBC->BCValue(mstime);
        nextBC=(NodalTempBC *)nextBC->GetNextObject();
    }
}

// Task 1b - get gradients in Vp * cp on particles
TransportTask *ConductionTask::GetGradients(double stepTime)
{
	int numnds,nds[maxShapeNodes];
    double fn[maxShapeNodes],xDeriv[maxShapeNodes],yDeriv[maxShapeNodes],zDeriv[maxShapeNodes];
	
	// Find gradients on the nonrigid particles
#pragma omp parallel for private(numnds,nds,fn,xDeriv,yDeriv,zDeriv)
    for(int p=0;p<nmpmsNR;p++)
	{	// find shape functions and derviatives
		const ElementBase *elref = theElements[mpm[p]->ElemID()];
		elref->GetShapeGradients(&numnds,fn,nds,mpm[p]->GetNcpos(),xDeriv,yDeriv,zDeriv,mpm[p]);
		
		// Find gradients from current temperatures
		mpm[p]->AddTemperatureGradient();			// zero gradient on the particle
		for(int i=1;i<=numnds;i++)
		{	Vector deriv = MakeVector(xDeriv[i],yDeriv[i],zDeriv[i]);
			mpm[p]->AddTemperatureGradient(ScaleVector(&deriv,nd[nds[i]]->gTemperature));
		}
	}
	
	return nextTask;
}

#pragma mark GRID FORCES EXTRAPOLATIONS

// find forces for conduction calculation (N-mm/sec = mJ/sec) (non-rigid particles only)
TransportTask *ConductionTask::AddForces(double *fcondBuffer,MPMBase *mptr,double sh,double dshdx,
										 double dshdy,double dshdz,TransportProperties *t)
{
	// internal force based on conduction tensor
	*fcondBuffer =  mptr->FCond(dshdx,dshdy,dshdz,t);
	
	// add source terms
	
	// if coupled to material dissipated energy, add and then zero dissipated energy
	if(energyCoupling)
	{	// V * q heat energy is mp (g) * specific energy (uJ/g) = uJ
		// To get mJ/sec, divide timestep (sec) and times 1e-3
		*fcondBuffer += sh*1.0e-3*mptr->mp*mptr->GetDispEnergy()/timestep;
	}
	
	// next task
	return nextTask;
}

// copy force from buffer to node
TransportTask *ConductionTask::AddForcesFromBuffer(NodalPoint *ndptr,double fcond)
{	ndptr->fcond += fcond;
	return nextTask;
}

// adjust forces at grid points with temperature BCs to have rates be correct
// to carry extrapolated temperatures (before impose BCs) to the correct
// one selected by grid based BC
TransportTask *ConductionTask::SetTransportForceBCs(double deltime)
{
    // Paste back noBC temperature
    int i;
    NodalTempBC *nextBC=firstTempBC;
    while(nextBC!=NULL)
        nextBC=nextBC->PasteNodalTemperature(nd[nextBC->GetNodeNum()]);
    
    // Set force to - mp Cp T(no BC)/timestep
	double mstime=1000.*(mtime+deltime);
    nextBC=firstTempBC;
    while(nextBC!=NULL)
	{   i=nextBC->GetNodeNum(mstime);
		if(i!=0) nd[i]->fcond = -nd[i]->gMpCp*nd[i]->gTemperature/deltime;
        nextBC=(NodalTempBC *)nextBC->GetNextObject();
	}
    
    // Now add each superposed concentration (* mp Cp) BC at incremented time
    nextBC=firstTempBC;
    while(nextBC!=NULL)
    {	i=nextBC->GetNodeNum(mstime);
		if(i!=0) nd[i]->fcond += nd[i]->gMpCp*nextBC->BCValue(mstime)/deltime;
        nextBC=(NodalTempBC *)nextBC->GetNextObject();
    }
	
	// --------- heat flux BCs -------------
	
	return nextTask;
}

#pragma mark UPDATE MOMENTA TASK

// get temperature rates on the nodes
TransportTask *ConductionTask::TransportRates(double deltime)
{
	// convert forces to temperature rates
	int i;
    for(i=1;i<=nnodes;i++)
	{   if(nd[i]->NumberNonrigidParticles()>0)
			nd[i]->fcond /= nd[i]->gMpCp;
	}
	return nextTask;
}
		
#pragma mark UPDATE PARTICLES TASK

// increment temperature rate on the particle
TransportTask *ConductionTask::IncrementTransportRate(NodalPoint *ndpt,double shape,double &rate) const
{	rate+=ndpt->fcond*shape;			// fcond are temperature rates from TransportRates()
	return nextTask;
}

// increment particle concentration (time is always timestep)
TransportTask *ConductionTask::MoveTransportValue(MPMBase *mptr,double deltime,double rate) const
{	mptr->pTemperature += deltime*rate;
	return nextTask;
}

#pragma mark UPDATE PARTICLE STRAIN TASK

// return increment transport rate
double ConductionTask::IncrementValueExtrap(NodalPoint *ndpt,double shape) const
{	return ndpt->gTemperature*shape;
}

// after extrapolated, find change this update on particle
double ConductionTask::GetDeltaValue(MPMBase *mptr,double pTempExtrap) const
{	double dTemperature = pTempExtrap-mptr->pPreviousTemperature;
	mptr->pPreviousTemperature = pTempExtrap;
	return dTemperature;
}

#pragma mark UPDATE STRAIN LAST TASK

// if needed for SZS or USAVG, update temperature on the grid (tempTime is always timestep)
TransportTask *ConductionTask::UpdateNodalValues(double tempTime)
{
	// add for each node
	int i;
    for(i=1;i<=nnodes;i++)
	{   if(nd[i]->NumberNonrigidParticles()>0)
			nd[i]->gTemperature += nd[i]->fcond*tempTime;
	}
	return nextTask;
}

#pragma mark CUSTOM METHODS

// If crack tip heating activated and there are cracks
// add heat of each crack as point sources for ndpt->fcond
void ConductionTask::AddCrackTipHeating(void)
{
	if(!crackTipHeating || firstCrack==NULL) return;
	CrackHeader *nextCrack=firstCrack;
	while(nextCrack!=NULL)
	{	nextCrack->CrackTipHeating();
		nextCrack=(CrackHeader *)nextCrack->GetNextObject();
	}
}

// Tell crack tip to heat itself when it propagates
void ConductionTask::StartCrackTipHeating(CrackSegment *crkTip,Vector &grow,double thickness)
{
	if(!crackTipHeating) return;
	double dist=sqrt(grow.x*grow.x+grow.y*grow.y);
	crkTip->StartCrackTipHeating(dist,thickness);
}



