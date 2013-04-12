/********************************************************************************
    ConductionTask.cpp
    NairnMPM
    
    Created by John Nairn on Fri Oct 15 2004
    Copyright (c) 2004 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Custom_Tasks/ConductionTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/MaterialBase.hpp" 
#include "Boundary_Conditions/NodalTempBC.hpp"
#include "Cracks/CrackHeader.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Cracks/CrackSegment.hpp"
#include "Materials/MicrostructureModel.hpp" //modiftf 
#include "Read_MPM/RPM.hpp" //modiftf #hardcodedheat

// global
bool ConductionTask::frictionalHeating=FALSE; //modiftf #frictionalheating
bool ConductionTask::hardCodedHeat=FALSE; //modiftf #hardcodedheat
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
{	double Cp=theMaterials[mptr->MatID()]->GetHeatCapacity(mptr);
	double rho;														
	if(rhoScaling)													//modiftf Vincent
		rho=(theMaterials[mptr->MatID()]->rho)/rhoScale;			//modiftf Vincent
	else															//modiftf Vincent
		rho=theMaterials[mptr->MatID()]->rho;						// original line
	double arg=mptr->volume*rho*Cp*shape;
	ndpt->gTemperature+=mptr->pTemperature*arg;
	ndpt->gRhoVCp+=arg;
	return nextTask;
}

// Task 1b - get grid temperatures and impose grid-based concentration BCs
void ConductionTask::GetValues(double stepTime)
{
	// convert to actual temperatires
	int i;
    for(i=1;i<=nnodes;i++)
	{   if(nd[i]->NumberNonrigidParticles()>0)
			nd[i]->gTemperature/=nd[i]->gRhoVCp;
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
	
	// impose flux boundary conditions (only allows zero for now)
}

// find forces for conduction calculation
TransportTask *ConductionTask::AddForces(NodalPoint *ndpt,MPMBase *mptr,double sh,double dshdx,double dshdy,double dshdz)
{
	// get conduction constant for this material (xx, yy, xy order)
	//Tensor *kten=theMaterials[mptr->MatID()]->GetkCondTensor();
	
	// internal force
	ndpt->fcond+=mptr->FCond(dshdx,dshdy,dshdz);
	//ndpt->fcond-=mptr->volume*((kten->xx*mptr->pTemp->DT.x + kten->xy*mptr->pTemp->DT.y)*dshdx
	//					+ (kten->xy*mptr->pTemp->DT.x + kten->yy*mptr->pTemp->DT.y)*dshdy);
	
	// add source terms
	
	// if coupled to material dissipated energy, add and then zero dissipated energy
	if(energyCoupling)
	{	// V * q = 1000 J/sec where J = 1.0e-9 * V (mm^3) * rho (g/cm^3) * specific energy
		// ...see data archiving for details on units
		double rho;
		if(rhoScaling)													// modiftf Vincent
			//rho=(theMaterials[mptr->MatID()]->rho)/rhoScale;		// modiftf Vincent
			rho=(theMaterials[mptr->MatID()]->rho)/1e6;
		else															// modiftf Vincent
			rho=theMaterials[mptr->MatID()]->rho;				// Original Line
		ndpt->fcond+=sh*1.0e-6*mptr->volume*rho*mptr->GetDispEnergy()/timestep;
	}
	
	//modiftf #frictionalheating
	if(frictionalHeating)
	{
		double rho=theMaterials[mptr->MatID()]->rho;
		
		//Solution based on particle mass:
		//ndpt->fcond+=sh*1.0e-6*mptr->volume*rho*ndpt->GetNodalFrictionEnergy()/timestep;

		//Solution based on material constant volume heat capacity:
		//ndpt->fcond+=sh*(1/(mptr->volume*rho*theMaterials[mptr->MatID()]->GetHeatCapacity(mptr)))*ndpt->GetNodalFrictionEnergy();
		
		//old
		//ndpt->fcond+=(sh*ndpt->GetNodalFrictionEnergy())/timestep;
	//debugtf 
	/*double datam = ndpt->GetNodalFrictionEnergy();
	if(datam>0)
		cout << "volume " << mptr->volume << endl; =6.25
		cout << "Node " << ndpt->num << " Friction Energy " << datam << endl;*/
	//debugtf modiftf 
		
		// based on actual needed values from John Nairn:
		// V * q * S, where q (J/s-cm^3), V(mm^3), S(dimensionless).
		//ndpt->fcond+=(sh*ndpt->GetNodalFrictionEnergy())/(timestep*1e3); //volumes cancel out with remaining 1e-3.
		//ndpt->fcond+=(ndpt->GetNodalFrictionEnergy()*1e2)/(timestep); //no shape function as from node already. volumes cancel out with remaining 1e2. (old)
		ndpt->fcond+=(ndpt->GetNodalFrictionEnergy()*1e2)/timestep; //no shape function as from node already. volumes cancel out with remaining 1e2. q is is J/s. (new)
		//store friction energy on material point?
		mptr->AddParticleFrictionEnergy(ndpt->GetNodalFrictionEnergy()); // actually nodal values
   
	}
	//modiftf #frictionalheating
	
	//modiftf #hardcodedheat
	if(hardCodedHeat)
	{
		// Variable power input
			double revpm = 900;						// Change this for various radial velocities
			double WW = (revpm*2*PI_CONSTANT)/60;
		
		// find tore from shear stress:
			Tensor *particleStress = mptr->GetStressTensor();
			//double Tore = particleStress->yy; //what is the shear stress term?
			//double Tore = 13.2e6;
			double rho=theMaterials[mptr->MatID()]->rho;
			double Tore = mptr->yieldC*rho/sqrt(3); // The shear yield stress from Microstructure Model(Pascals)
		
		// find position of particle relative to tool centre
			double X=mptr->pos.x-rotator->getxCentre();		// to store relative x position
			double Y=mptr->pos.y-rotator->getyCentre();		// to store relative y position
			//double Z=mptr->pos.z;							// to store relative z position
			double Z=mptr->pos.z-rotator->getzDepth();
			
			// for Convective Flux
			double XX = mptr->pos.x;
			double YY = mptr->pos.y;
			
			// For Radial Input Power
			double r2 = (pow(Y,2)+pow(X,2));
			double r1 = sqrt(r2)*(1.0e-1)/(0.1);  // account for radius stored in mm->cm and particle spacing (cm).

			// for varibale power use below:
			double inpower = Tore*WW*1e-6*r1;
			
			
			
			// 
			if(Z>=-0.5&&Z<3.5)	 // pin heat flux
			{	if((r2)<=(pow(2.5,2))) // 0.5 + outside radius of pin EDIT FOR PIN RADIUS
				{	ndpt->fcond+=mptr->volume*inpower*sh;}
				//cout << "1. pos.z: " << mptr->pos.z << " tool z depth: " << rotator->getzDepth() << endl;}
			}
			
			else if(Z>=3.5)	// Just under Shoulder begins EDIT FOR HEIGHT OF PIN
			{	if((r2)<=(pow(6.5,2))) // Exact radius of Shoulder EDIT FOR SHOULDER RADIUS
				{	ndpt->fcond+=mptr->volume*inpower*sh;}
				//cout << "2. pos.z: " << mptr->pos.z << " tool z depth: " << rotator->getzDepth() << endl;}
			}
		
			
			// Convective Flux at sides and top -> 15W/m2 = 15e-4W/cm2 -> 15e-3W/cm^3. As this is a coefficient * temperature difference, where 22 is room temp.
			
			// version John Mikhail (JM)			
			//if(XX<1.5||XX>200.5||YY<0.5||YY>49.5||Z>7.5)
			//{	ndpt->fcond-=mptr->volume*0.015*sh*(mptr->pTemperature-24);	// V * q * S, where q (J/s-cm^3), V(mm^3), S(dimensionless).
				//cout << mptr->pos.x << " " << mptr->pos.y << " " << mptr->pos.z << endl;
			//}
			
				
			// Backing Plate
			// version John Mikhail (JM) and Sri Lathabai (SL)
			//if(Z<2.5) //<2
			//{	ndpt->fcond-=mptr->volume*0.5*sh*(mptr->pTemperature-24);	// V * q * S, where q (J/s-cm^3), V(mm^3), S(dimensionless).
				//cout << mptr->pos.x << " " << mptr->pos.y << " " << mptr->pos.z << endl;
			//}	
	
	}
	// modiftf #hardcodedheat
	
	// if on boundary, get boundary force
	
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
    
    // Set force to - rho V Cp T(no BC)/timestep
	double mstime=1000.*(mtime+deltime);
    nextBC=firstTempBC;
    while(nextBC!=NULL)
	{   i=nextBC->GetNodeNum(mstime);
		if(i!=0) nd[i]->fcond=-nd[i]->gRhoVCp*nd[i]->gTemperature/deltime;
        nextBC=(NodalTempBC *)nextBC->GetNextObject();
	}
    
    // Now add each superposed concentration (* rho V Cp) BC at incremented time
    nextBC=firstTempBC;
    while(nextBC!=NULL)
    {	i=nextBC->GetNodeNum(mstime);
		if(i!=0) nd[i]->fcond+=nd[i]->gRhoVCp*nextBC->BCValue(mstime)/deltime;
        nextBC=(NodalTempBC *)nextBC->GetNextObject();
    }
	
	return nextTask;
}

// get temperature rates on the nodes
TransportTask *ConductionTask::TransportRates(double deltime)
{
	// convert forces to temperature rates
	int i;
    for(i=1;i<=nnodes;i++)
	{   if(nd[i]->NumberNonrigidParticles()>0)
			nd[i]->fcond/=nd[i]->gRhoVCp;
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
			nd[i]->gTemperature+=nd[i]->fcond*tempTime;
	}
	return nextTask;
}

// increment transport rate
TransportTask *ConductionTask::IncrementValueExtrap(NodalPoint *ndpt,double shape)
{	pValueExtrap+=ndpt->gTemperature*shape;
	return nextTask;
}

// after extrapolated, find change this update on particle
TransportTask *ConductionTask::GetDeltaValue(MPMBase *mptr)
{	dTemperature=pValueExtrap-mptr->pPreviousTemperature;
	mptr->pPreviousTemperature=pValueExtrap;
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



