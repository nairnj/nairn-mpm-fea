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
		Divide gTemperature by gMpCp and impose grid T BCs (GetValues())
		Find grad T on particle (GetGradients())
	Grid Forces Task
		Extrapolate conductivity force to fcond (AddForces())
			(include heat sources, and energy coupling)
		Add crack tip heating to fcond (AddCrackTipHeating())
		Finish grid BCs and impose flux BCs in fcond (SetTransportForceBCs())
	Update Momenta Task
		Divide fcond by gMpCp to get temperature rates (TransportRates())
	Update Particles Task
		Each particle: zero rate, extrapolate from nodes to particle, then
			update particle (ZeroTransportRate(), IncrementTransportRate(),
			MoveTransportRate()).
	Update strains last (if used)
		Update gTemperature on nodes (UpdateNodalValues())
	Update strains on particles (both places)
		Extrapolate grid temperature to particle using pValueExtra (IncrementValueExtrap())
		Find dTemperature = extrapolated value minus previous extrapolated value and then
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

// global
bool ConductionTask::active=FALSE;
bool ConductionTask::crackTipHeating=FALSE;
bool ConductionTask::energyCoupling=FALSE;
double ConductionTask::dTemperature=0.;
ConductionTask *conduction=NULL;

#pragma mark INITIALIZE

// Constructors
ConductionTask::ConductionTask()
{	// allocate diffusion data on each particle
    int p;
	for(p=0;p<nmpms;p++)
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

// adjust time for given cell size if needed
TransportTask *ConductionTask::TransportTimeStep(int matid,double dcell,double *tmin)
{	double diffCon=theMaterials[matid]->MaximumDiffusivity();
	double tst=(dcell*dcell)/(4.*diffCon);						// factor 2 shorter than minimum
	if(tst<*tmin) *tmin=tst;
	return nextTask;
}

#pragma mark TASK EXTRAPOLATION METHODS

// Task 1 Extrapolation of concentration to the grid
TransportTask *ConductionTask::Task1Extrapolation(NodalPoint *ndpt,MPMBase *mptr,double shape)
{	double Cp=theMaterials[mptr->MatID()]->GetHeatCapacity(mptr);   // mJ/(g-K)
	double arg = mptr->mp*Cp*shape;                                 // mJ/K
	ndpt->gTemperature+=mptr->pTemperature*arg;                     // mJ
	ndpt->gMpCp+=arg;                                               // mJ/K
	return nextTask;
}

// Task 1b - get grid temperatures and impose grid-based concentration BCs
void ConductionTask::GetValues(double stepTime)
{
	// convert to actual temperatures
	int i;
    for(i=1;i<=nnodes;i++)
	{   if(nd[i]->NumberNonrigidParticles()>0)
			nd[i]->gTemperature /= nd[i]->gMpCp;
	}

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
void ConductionTask::GetGradients(double stepTime)
{
    int i,p,iel;
    double fn[MaxShapeNds],xDeriv[MaxShapeNds],yDeriv[MaxShapeNds],zDeriv[MaxShapeNds];
	int numnds,nds[MaxShapeNds];
	
	// Find gradients on the particles
    for(p=0;p<nmpms;p++)
	{	if(theMaterials[mpm[p]->MatID()]->Rigid()) continue;
	
		// find shape functions and derviatives
		iel=mpm[p]->ElemID();
		theElements[iel]->GetShapeGradients(&numnds,fn,nds,mpm[p]->GetNcpos(),xDeriv,yDeriv,zDeriv,mpm[p]);
		
		// Find gradients from current temperatures
		mpm[p]->AddTemperatureGradient();			// zero gradient on the particle
		for(i=1;i<=numnds;i++)
		{	Vector deriv=MakeVector(xDeriv[i],yDeriv[i],zDeriv[i]);
			mpm[p]->AddTemperatureGradient(ScaleVector(&deriv,nd[nds[i]]->gTemperature));
		}
	}
}

// find forces for conduction calculation (N-mm/sec = mJ/sec)
TransportTask *ConductionTask::AddForces(NodalPoint *ndpt,MPMBase *mptr,double sh,double dshdx,double dshdy,double dshdz)
{
	// internal force based on conduction tensor
	ndpt->fcond += mptr->FCond(dshdx,dshdy,dshdz);
	
	// add source terms
	
	// if coupled to material dissipated energy, add and then zero dissipated energy
	if(energyCoupling)
	{	// V * q heat energy is mp (g) * specific energy (uJ/g) = uJ
		// To get mJ/sec, divide timestep (sec) and times 1e-3
		ndpt->fcond += sh*1.0e-3*mptr->mp*mptr->GetDispEnergy()/timestep;
	}
	
	// next task
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
		
// increment temperature rate on the particle
TransportTask *ConductionTask::IncrementTransportRate(NodalPoint *ndpt,double shape)
{	rate+=ndpt->fcond*shape;			// fcond are temperature rates from TransportRates()
	return nextTask;
}

// increment particle concentration (time is always timestep)
TransportTask *ConductionTask::MoveTransportValue(MPMBase *mptr,double deltime)
{	mptr->pTemperature+=deltime*rate;
	return nextTask;
}

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

// increment transport rate
TransportTask *ConductionTask::IncrementValueExtrap(NodalPoint *ndpt,double shape)
{	pValueExtrap += ndpt->gTemperature*shape;
	return nextTask;
}

// after extrapolated, find change this update on particle
TransportTask *ConductionTask::GetDeltaValue(MPMBase *mptr)
{	dTemperature = pValueExtrap-mptr->pPreviousTemperature;
	mptr->pPreviousTemperature = pValueExtrap;
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



