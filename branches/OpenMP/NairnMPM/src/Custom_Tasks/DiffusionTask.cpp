/********************************************************************************
    DiffusionTask.cpp
    NairnMPM
    
    Created by John Nairn on Mon Mar 08 2004
    Copyright (c) 2003 John A. Nairn, All rights reserved.
 
    Diffusion calculations
   -------------------------
    See comments in ConductionTask.cpp but:
    Change gTemperature to gConcentration, gMpCp to gVolume, and fcond to fdiff
    Update Particles Task
        cut off particle potential to range 0 to 1
    Chemical potential (or concentration potential 0 to 1)
        Internally all calculations in terms or potential (0 to 1)
        Output concentration and concentration gradient scaled to
            csat to get weight fraction and weight fraction gradient
        Flux set as mass per area per sec. See documentation for conversions.
********************************************************************************/

#include "Custom_Tasks/DiffusionTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "Boundary_Conditions/NodalConcBC.hpp"
#include "Boundary_Conditions/MatPtFluxBC.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Nodes/NodalPoint.hpp"

// NEWINCLUDE
#include "Exceptions/CommonException.hpp"

// global
bool DiffusionTask::active=FALSE;
double DiffusionTask::reference = 0.;				// zero-strain concentration
DiffusionTask *diffusion=NULL;

#pragma mark INITIALIZE

// Constructors
DiffusionTask::DiffusionTask()
{	// allocate diffusion data on each particle
    // done before know number of nonrigid, so do on all
	for(int p=0;p<nmpms;p++)
		mpm[p]->AllocateDiffusion();
}

// Return name of this task
const char *DiffusionTask::TaskName(void) { return "diffusion analysis"; }

#pragma mark STANDARD METHODS

// diffsion analysis settings
TransportTask *DiffusionTask::TransportOutput(void)
{	char mline[200];
	TransportTask::TransportOutput();
	sprintf(mline,"   Reference concentration =%8.4lf",reference);
	cout << mline << endl;
	return nextTask;
}

// adjust time for given cell size if needed
TransportTask *DiffusionTask::TransportTimeStep(int matid,double dcell,double *tmin)
{	double diffCon=theMaterials[matid]->MaximumDiffusion();
	double tst=(dcell*dcell)/(4.*diffCon);						// factor 2 shorter than minimum
	if(tst<*tmin) *tmin=tst;
	return nextTask;
}

#pragma mark TASK EXTRAPOLATION METHODS

// Task 1 Extrapolation of concentration to the grid. Concentration is acutally a chemical
// potential from 0 to 1 where 1 means concentration is equal to that materials saturation
// concentration. (units are mm^3)
TransportTask *DiffusionTask::Task1Extrapolation(NodalPoint *ndpt,MPMBase *mptr,double shape)
{   double Vp = mptr->GetVolume(DEFORMED_VOLUME);
	ndpt->gConcentration += mptr->pConcentration*Vp*shape;
	ndpt->gVolume += Vp*shape;
	return nextTask;
}

// Task 1 reduction of ghost node to real node for concentration on the grid
TransportTask *DiffusionTask::Task1Reduction(NodalPoint *ndpt,NodalPoint *ghostNdpt)
{	ndpt->gConcentration += ghostNdpt->gConcentration;                    // mJ
	ndpt->gVolume += ghostNdpt->gVolume;                                  // mJ/K
	return nextTask;
}

// Get grid concentrations and impose grid-based concentration BCs
TransportTask *DiffusionTask::GetNodalValue(NodalPoint *ndptr)
{	if(ndptr->NumberNonrigidParticles()>0)
		ndptr->gConcentration /= ndptr->gVolume;
	return nextTask;
}

// Get grid concentrations and impose grid-based concentration BCs
void DiffusionTask::ImposeValueBCs(double stepTime)
{
	int i;
	
	// Copy no-BC concentration
    NodalConcBC *nextBC=firstConcBC;
    while(nextBC!=NULL)
		nextBC=nextBC->CopyNodalConcentration(nd[nextBC->GetNodeNum()]);
	
    // Set active ones to zero
	double mstime=1000.*stepTime;
    nextBC=firstConcBC;
    while(nextBC!=NULL)
	{   i=nextBC->GetNodeNum(mstime);
	    if(i!=0) nd[i]->gConcentration = 0.;
        nextBC=(NodalConcBC *)nextBC->GetNextObject();
    }

    // Now add all concentrations to nodes with active concentration BCs
    nextBC=firstConcBC;
    while(nextBC!=NULL)
	{   i=nextBC->GetNodeNum(mstime);
	    if(i!=0) nd[i]->gConcentration += nextBC->BCValue(mstime);
        nextBC=(NodalConcBC *)nextBC->GetNextObject();
    }
	
	// verify all set BCs are between 0 and 1
    nextBC=firstConcBC;
    while(nextBC!=NULL)
	{   i=nextBC->GetNodeNum(mstime);
		if(i!=0)
		{	if(nd[i]->gConcentration<0)
				nd[i]->gConcentration = 0.;
			else if(nd[i]->gConcentration>1.)
				nd[i]->gConcentration = 1.;
		}
        nextBC=(NodalConcBC *)nextBC->GetNextObject();
    }
}

// Get gradients in Vp * cp on particles
TransportTask *DiffusionTask::GetGradients(double stepTime)
{
    CommonException *transErr = NULL;
	int nds[maxShapeNodes];
    double fn[maxShapeNodes],xDeriv[maxShapeNodes],yDeriv[maxShapeNodes],zDeriv[maxShapeNodes];
	
	// Find gradients on the nonrigid particles
#pragma omp parallel for private(nds,fn,xDeriv,yDeriv,zDeriv)
    for(int p=0;p<nmpmsNR;p++)
	{	try
        {   // find shape functions and derviatives
            MPMBase *mptr = mpm[p];
            const ElementBase *elref = theElements[mptr->ElemID()];
            int i,numnds;
            elref->GetShapeGradients(&numnds,fn,nds,mptr->GetNcpos(),xDeriv,yDeriv,zDeriv,mptr);
            
            // Find gradients from current concentrations
            mptr->AddConcentrationGradient();			// zero gradient on the particle
            for(i=1;i<=numnds;i++)
            {	Vector deriv=MakeVector(xDeriv[i],yDeriv[i],zDeriv[i]);
                mptr->AddConcentrationGradient(ScaleVector(&deriv,nd[nds[i]]->gConcentration));
            }
        }
        catch(CommonException err)
        {   if(transErr!=NULL)
            {
#pragma omp critical
                transErr = new CommonException(err);
            }
        }
	}
    
    // throw any errors
    if(transErr!=NULL) throw *transErr;
	
	return nextTask;
}

// find forces for diffusion calculation (mm^3/sec) (non-rigid particles only)
TransportTask *DiffusionTask::AddForces(NodalPoint *ndptr,MPMBase *mptr,double sh,double dshdx,
										double dshdy,double dshdz,TransportProperties *t)
{
	// internal force
	ndptr->fdiff += mptr->FDiff(dshdx,dshdy,dshdz,t);
	
	// add source terms (should be potential per sec, if c units per second, divide by concSaturation)
	
	// return next task
	return nextTask;
}

// copy diffusion forces from ghost node to real node
TransportTask *DiffusionTask::CopyForces(NodalPoint *ndptr,NodalPoint *ghostNdptr)
{	ndptr->fdiff += ghostNdptr->fdiff;
	return nextTask;
}

// adjust forces at grid points with concentration BCs to have rates be correct
// to carry extrapolated concentrations (before impose BCs) to the correct
// one selected by grid based BC
TransportTask *DiffusionTask::SetTransportForceBCs(double deltime)
{
    int i;
    NodalConcBC *nextBC=firstConcBC;
    
	// --------- consistent forces for grid concentration BCs ------------
	
    // Paste back noBC concentration
    while(nextBC!=NULL)
        nextBC = nextBC->PasteNodalConcentration(nd[nextBC->GetNodeNum()]);
    
    // Set force to - VC(no BC)/timestep
	double mstime=1000.*(mtime+deltime);
    nextBC=firstConcBC;
    while(nextBC!=NULL)
    {	i = nextBC->GetNodeNum(mstime);
		if(i!=0) nd[i]->fdiff = -nd[i]->gVolume*nd[i]->gConcentration/deltime;
		nextBC = (NodalConcBC *)nextBC->GetNextObject();
    }

    // Now add each superposed concentration (* volume) BC at incremented time
    nextBC=firstConcBC;
    while(nextBC!=NULL)
    {	i = nextBC->GetNodeNum(mstime);
		if(i!=0) nd[i]->fdiff += nd[i]->gVolume*nextBC->BCValue(mstime)/deltime;
        nextBC = (NodalConcBC *)nextBC->GetNextObject();
    }
	
	// --------- concentration flux BCs -------------
	
	MatPtFluxBC *nextFlux=firstFluxPt;
    while(nextFlux!=NULL)
    	nextFlux = nextFlux->AddMPFlux(deltime);

	return nextTask;
}

// find concentration rates on the nodes, but since these are potentials they are limited
// to 0 to 1 or limited to 0 < gConcentration + fdiff*deltime < 1
TransportTask *DiffusionTask::TransportRates(NodalPoint *ndptr,double deltime)
{
	if(ndptr->NumberNonrigidParticles()>0)
	{	ndptr->fdiff /= ndptr->gVolume;
		/*
		double concTest = ndptr->gConcentration + ndptr->fdiff*deltime;
		if(concTest<0.)
			ndptr->fdiff = -ndptr->gConcentration/deltime;          // will evolve to 0
		else if(concTest>1.)
			ndptr->fdiff = (1.-ndptr->gConcentration)/deltime;      // will evolve to 1
		*/
	}
	return nextTask;
}
		
// increment concentration rate on the particle
TransportTask *DiffusionTask::IncrementTransportRate(NodalPoint *ndpt,double shape,double &rate) const
{	rate += ndpt->fdiff*shape;			// fdiff are concentration rates from TransportRates()
	return nextTask;
}

// increment particle concentration
TransportTask *DiffusionTask::MoveTransportValue(MPMBase *mptr,double deltime,double rate) const
{	mptr->pConcentration += deltime*rate;
    if(mptr->pConcentration<0.)
        mptr->pConcentration = 0.;
    else if(mptr->pConcentration>1.)
        mptr->pConcentration = 1.;
	return nextTask;
}

// if needed for SZS or USAVG, update concentration on the grid (concTime is always timestep)
TransportTask *DiffusionTask::UpdateNodalValues(double concTime)
{
	// add for each node
	int i;
    for(i=1;i<=nnodes;i++)
	{   if(nd[i]->NumberNonrigidParticles()>0)
			nd[i]->gConcentration += nd[i]->fdiff*concTime;
	}
	return nextTask;
}

// increment transport rate
double DiffusionTask::IncrementValueExtrap(NodalPoint *ndpt,double shape) const
{	return ndpt->gConcentration*shape;
}

// after extrapolated, find change this update on particle and reset particle
// property to this grid extrapolated value
double DiffusionTask::GetDeltaValue(MPMBase *mptr,double pConcExtrap) const
{	double dConcentration = pConcExtrap-mptr->pPreviousConcentration;
	mptr->pPreviousConcentration = pConcExtrap;
	return dConcentration;
}

		
