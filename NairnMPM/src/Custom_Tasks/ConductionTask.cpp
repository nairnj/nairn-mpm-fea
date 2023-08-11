/*********************************************************************************************
    ConductionTask.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Fri Oct 15 2004
    Copyright (c) 2004 John A. Nairn, All rights reserved.
 
	Conduction calculations
   -------------------------
	See TransportTask.cpp comments with
		gTValue, gVCT, and gQ in gCond
		gV is used only in revised heat method
 
	Unique things
    --------------
	Grid Forces Task
		Add crack tip heating to gCond.gQ (AddCrackTipHeating())
	Update Particles Task
		Add adiabatic heating to particles
*********************************************************************************************/

#include "stdafx.h"
#include "Custom_Tasks/ConductionTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "Boundary_Conditions/NodalTempBC.hpp"
#include "Boundary_Conditions/MatPtHeatFluxBC.hpp"
#include "Cracks/CrackHeader.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Cracks/CrackSegment.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Materials/RigidMaterial.hpp"
#include "Exceptions/CommonException.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "System/UnitsController.hpp"

// flag to activate
bool ConductionTask::active=false;
bool ConductionTask::activeRamp=false;

// flags only allowed when conduction is active
bool ConductionTask::crackTipHeating=false;
bool ConductionTask::crackContactHeating=false;
bool ConductionTask::matContactHeating=false;

// these don't require conduction to be active
bool ConductionTask::adiabatic=false;

// possible one conduction task
ConductionTask *conduction=NULL;

#pragma mark STANDARD METHODS

// Return name of this task
const char *ConductionTask::TaskName(void) { return "conduction calculations"; }

// called once at start of MPM analysis and after preliminary calcs are done
TransportTask *ConductionTask::Initialize(void)
{
	// gradient addresses if needed
	crackGradT = -1;
	materialGradT = -1;
	
	// print task details
	cout << "Coupled " << TaskName() << endl;
	
	// time step
	char fline[256];
	size_t fsize=256;
	snprintf(fline,fsize,"   Conduction time step maximum (%s): %.7e",UnitsController::Label(ALTTIME_UNITS),transportTimeStep*UnitsController::Scaling(1.e3));
	cout << fline << endl;
	cout << "   Time step factor: " << fmobj->GetTransCFLCondition() << endl;
	
	// features
	if(crackTipHeating)
		cout << "   Crack tip heating activated" << endl;
	if(crackContactHeating)
		cout << "   Crack contact frictional heating activated" << endl;
	if(matContactHeating)
		cout << "   Material contact frictional heating activated" << endl;
	
	// allocate conduction data on each particle (reservoir too)
    // done before know number of nonrigid, so do on all
	int size = 3;
	if(crackGradT>0) size+=3;
	if(materialGradT>0) size+=3;
	for(int p=0;p<nmpms;p++)
	{	mpm[p]->pTemp = new double[size];
		for(int i=0;i<size;i++) mpm[p]->pTemp[i] = 0.;
	}
	
	return nextTask;
}

#pragma mark MASS AND MOMENTUM EXTRAPOLATIONS

// Task 1 Extrapolation of temperature to the grid
// Only called for non-rigid materials
TransportTask *ConductionTask::Task1Extrapolation(NodalPoint *ndpt,MPMBase *mptr,double shape,short vfld,int matfld)
{
	double Cv = theMaterials[mptr->MatID()]->GetHeatCapacity(mptr);		// nJ/(g-K) using Cv is correct
	double CTShape = mptr->mp*Cv*shape;
	TransportField *gTrans = GetTransportFieldPtr(ndpt);
	double mTpTValueShape = GetParticleValue(mptr)*CTShape;
	gTrans->gTValue += mTpTValueShape;
	gTrans->gVCT += CTShape;
	Task1ContactExtrapolation(ndpt,vfld,matfld,mTpTValueShape,CTShape);
	return nextTask;
}

// Get Vp * CTp
double ConductionTask::GetVpCTp(MPMBase *mptr)
{	double Cv = theMaterials[mptr->MatID()]->GetHeatCapacity(mptr);		// nJ/(g-K) using Cv is correct
	return mptr->mp*Cv;
}

// Zero gradients on the particles
void ConductionTask::ZeroTransportGradients(MPMBase *mptr)
{	// zero gradient for global T on the particle
	mptr->AddTemperatureGradient(GRAD_GLOBAL);
    ZeroTransportContactGradients(mptr);
}

// Add gradients on the particles
void ConductionTask::AddTransportGradients(MPMBase *mptr,Vector *deriv,NodalPoint *ndptr,short vfld)
{   Vector TGip;
	mptr->AddTemperatureGradient(GRAD_GLOBAL,CopyScaleVector(&TGip,deriv,ndptr->gCond.gTValue));
	AddTransportContactGradients(mptr,deriv,ndptr,vfld);
}

#pragma mark GRID FORCES EXTRAPOLATIONS

// find forces for conduction calculation (N-mm/sec = mJ/sec) (non-rigid particles only)
TransportTask *ConductionTask::AddForces(NodalPoint *ndptr,MPMBase *mptr,double sh,double dshdx,
										 double dshdy,double dshdz,TransportProperties *t,short vfld,int matfld)
{	// internal force based on conduction tensor
	ndptr->gCond.gQ += mptr->FCond(GRAD_GLOBAL,dshdx,dshdy,dshdz,t);
	
	// add source terms
	
	// contact forces
	AddContactForces(ndptr,mptr,sh,dshdx,dshdy,dshdz,t,vfld,matfld);
	
	return nextTask;
}

#ifndef TRANSPORT_FMPM
// adjust forces at grid points with temperature BCs to have rates be correct
// to carry extrapolated temperatures (before impose BCs) to the correct
// one selected by grid based BC
// only used when TRANSPORT_FMPM is not define; that method uses new methods
TransportTask *ConductionTask::SetTransportForceAndFluxBCs(double deltime)
{
    // Paste back noBC temperature
    int i;
    NodalTempBC *nextBC=firstTempBC;
    while(nextBC!=NULL)
    {   i=nextBC->GetNodeNum(mtime);
		if(i!=0)
		{   nextBC->PasteNodalValue(nd[i],GetTransportFieldPtr(nd[i]));
			nextBC->InitQReaction();
			double qflow = -nd[i]->gCond.gQ;
			nextBC->SuperposeQReaction(qflow);
			nd[i]->gCond.gQ = 0.;
		}
        nextBC=(NodalTempBC *)nextBC->GetNextObject();
	}

    // Set force to - T(no BC)/timestep (only once per node)
    nextBC=firstTempBC;
    while(nextBC!=NULL)
	{   i=nextBC->GetNodeNum(mtime);
		if(i!=0)
		{	// but only once per node in case more than one Temperature BC on the node
			if(nd[i]->gCond.gQ==0.)
			{	// Power (energy/time)
				double qflow = -nd[i]->gCond.gVCT*nd[i]->gCond.gTValue/deltime;
				nd[i]->gCond.gQ = qflow;
				// for global archive of boundary heat
				nextBC->SuperposeQReaction(qflow);
			}
		}
        nextBC=(NodalTempBC *)nextBC->GetNextObject();
	}
    
    // Now add each superposed temperature BC at incremented time
	// Can superpose T, but one should be absolute T and others as T increments
    nextBC=firstTempBC;
    while(nextBC!=NULL)
    {	i=nextBC->GetNodeNum(mtime);
		if(i!=0)
		{	// Power (energy/time)
			double qflow = nd[i]->gCond.gVCT*nextBC->BCValue(mtime)/deltime;
			nd[i]->gCond.gQ += qflow;
			// for global archive of boundary flow
			nextBC->SuperposeQReaction(qflow);
		}
        nextBC=(NodalTempBC *)nextBC->GetNextObject();
    }
	
	// --------- heat flux BCs -------------
	MatPtLoadBC *nextFlux = firstHeatFluxPt;
    while(nextFlux!=NULL)
    	nextFlux = nextFlux->AddMPFluxBC(mtime);
	
	// next task
	return nextTask;
}
#endif

#pragma mark UPDATE PARTICLES TASKS

#pragma mark CUSTOM METHODS

// If crack tip heating activated and there are cracks
// add heat of each crack as point sources for ndpt->gCond.gQ
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

#pragma mark ACCESSORS

// adjust time for given cell size if needed
TransportTask *ConductionTask::TransportTimeStepFactor(int matid,double *diffCon)
{	*diffCon = theMaterials[matid]->MaximumDiffusivity();
    return nextTask;
}

// return point on node to transport field
TransportField *ConductionTask::GetTransportFieldPtr(NodalPoint *ndpt) const { return &(ndpt->gCond); }

// return first boundary condition
NodalValueBC *ConductionTask::GetFirstBCPtr(void) const { return firstTempBC; }
MatPtLoadBC *ConductionTask::GetFirstFluxBCPtr(void) const { return firstHeatFluxPt; }

// particle values
double ConductionTask::GetParticleValue(MPMBase *mptr) const { return mptr->pTemperature; }
double *ConductionTask::GetParticleValuePtr(MPMBase *mptr) const { return &(mptr->pTemperature); }
double ConductionTask::GetPrevParticleValue(MPMBase *mptr) const { return mptr->pPreviousTemperature; }
double *ConductionTask::GetPrevParticleValuePtr(MPMBase *mptr) const { return &(mptr->pPreviousTemperature); }

#pragma mark CLASS METHODS

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
    if(adiabatic)
        cout << " and adiabatic" << endl;
    else
        cout << " and isothermal" << endl;
}

// is the system isolated?
// update this if new thermal BCs are added
bool ConductionTask::IsSystemIsolated(void)
{
    // if has ramp, then is it not isolated
	if(activeRamp) return FALSE;
    
    // if no conduction then is isolated
    if(!active) return TRUE;
    
    // if conduction is active, still isolated if no thermal BC
    if(firstTempBC!=NULL) return FALSE;
    
    // check if any active rigid particles set temperature
	if(RigidMaterial::someSetTemperature) return FALSE;
    
    // must be isolated
    return TRUE;
}
