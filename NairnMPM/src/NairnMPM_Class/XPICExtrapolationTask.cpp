/********************************************************************************
	XPICExtrapolationTask.cpp
	nairn-mpm-fea

	Created by John Nairn on June 30, 2016
	Copyright (c) 2016 John A. Nairn, All rights reserved.
 
	The task only activated with XPIC is on.
 
	The tasks are:
	--------------
	* Find vstar on all nodes (needed when update particles)
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/XPICExtrapolationTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "Exceptions/CommonException.hpp"
#include "Patches/GridPatch.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Cracks/CrackNode.hpp"
#include "Nodes/MaterialContactNode.hpp"

// class globals
XPICExtrapolationTask *XPICMechanicsTask=NULL;

#pragma mark CONSTRUCTORS

XPICExtrapolationTask::XPICExtrapolationTask(const char *name) : MPMTask(name)
{
	// It always does FMPM(2) at least
	dynamicKmax = 2;
}

#pragma mark REQUIRED METHODS

// XPIC and FMPM extrapolation to get modifed grid velocities or transport values
// xpicOption is XPIC_STRAIN_UPDATE or XPIC_PARTICLE_UPDATE for two ways this method
//		is called. First only occurs in USF, USAVG+, USAVG-, and USF. It is 0
//		when called for transpoirt calculations
// WARNING: this entry point is called by both mechanics and transport. Any features
//		added here about mechanics should be overridden in TO version
// throws CommonException()
bool XPICExtrapolationTask::Execute(int xpicOption)
{
	// get order for this version of XPIC and skip if no calculations needed
	int m = GetXPICOrder();
	if(m==0) return true;
	
	CommonException *xpicErr = NULL;
#ifdef CONST_ARRAYS
	double fn[MAX_SHAPE_NODES];
	int ndsArray[MAX_SHAPE_NODES];
#else
	double fn[maxShapeNodes];
	int ndsArray[maxShapeNodes];
#endif
	
	// For Mechanics:
	//   Set v* = v^{L+} for FMPM or v^{L} for XPIC (in end v(k)=v*)
	//       (for XPIC add a*dt to v* so result will be v(k)+a*dt)
	//   Set v*(prev) = v* on real nodes
	//   Set v*(next) = 0 on real and ghost nodes
	// For Transport set corresponding terms for transport value
	//   Always FMPM for transport
#pragma omp parallel for
	for(int i=1;i<=*nda;i++)
		InitializeXPICData(nd[nda[i]],timestep,xpicOption);
	
	// In parallel, initialize patches too
	int totalPatches = fmobj->GetTotalNumberOfPatches();
	if(totalPatches>1)
	{	for(int i=0;i<totalPatches;i++)
			InitializeXPICData(patches[i],xpicOption);
	}
	
	// iterate for k from 2 to find v_k*
	//   (note that loop is skipped for order 1)
	for(int k=2;k<=m;k++)
	{
#pragma omp parallel private(fn,ndsArray)
		{
			// thread for patch pn
			int pn = GetPatchNumber();
			
			try
			{	// Loop over non-rigid particles
				MPMBase *mpmptr = patches[pn]->GetFirstBlockPointer(FIRST_NONRIGID);
				while(mpmptr!=NULL)
				{	// get shape functions and mat field
					const MaterialBase *matID = theMaterials[mpmptr->MatID()];			// material object for this particle
					int matfld = matID->GetField();										// material velocity field
					
					// get nodes and shape function for material point p
					const ElementBase *elref = theElements[mpmptr->ElemID()];
					int *nds = ndsArray;
					elref->GetShapeFunctions(fn,&nds,mpmptr);
					
					// double loop over nodes
					// Subtract (S^+S)v*(prev) from v*(next)
					// With velocity gradient, subtract (S^+S + S^{L+}S^L)v*(prev) from v*(next)
					XPICDoubleLoop(mpmptr,matfld,nds,fn,pn);
					
					// next material point
					mpmptr = (MPMBase *)mpmptr->GetNextObject();
				}
			}
			catch(CommonException& err)
			{   if(xpicErr==NULL)
				{
#pragma omp critical (error)
					xpicErr = new CommonException(err);
				}
			}
			catch(std::bad_alloc&)
			{   if(xpicErr==NULL)
				{
#pragma omp critical (error)
					xpicErr = new CommonException("Memory error","XPICExtrapolationTask::Execute");
				}
			}
			catch(...)
			{   if(xpicErr==NULL)
				{
#pragma omp critical (error)
					xpicErr = new CommonException("Unexpected error","XPICExtrapolationTask::Execute");
				}
			}
		}
		
		// back to serial

		// throw now - only possible error if too many CPDI nodes in 3D
		if(xpicErr!=NULL) throw *xpicErr;
		
		// reduction of ghost node velocities to real nodes
		// add ghost to real nodes and set vnext on ghost to zero for next time step.
		if(totalPatches>1)
		{	for(int pn=0;pn<totalPatches;pn++)
				ReduceXPICData(patches[pn]);
		}
		
		// Set Delta v(k) as vprev = vprev - vnext
		// Set BC nodes to zero and do incremental material contact
		GetDeltaV(timestep,xpicOption);

		// Increment v(k) += Delta v(k) (in code v* += vprev)
		// Zero vnext on real nodes for next pass through the loop
#pragma omp parallel for
		for(int i=1;i<=*nda;i++)
			UpdateXStar(nd[nda[i]],timestep);

	}

    return true;
}

// Get order to velocity XPIC calculations
int XPICExtrapolationTask::GetXPICOrder(void)
{	// The XPIC order (always 2 or greater when this task is activated)
	int m = bodyFrc.GetXPICOrder();
	
	// Skip if FLIP or PIC (=XPIC(1))
	return m>1 ? m : 0 ;
}

// Initialize real node vStar, vStarNext, and vStarPrev
void XPICExtrapolationTask::InitializeXPICData(NodalPoint *ndptr,double timestep,int xpicOption)
{	ndptr->XPICSupport(INITIALIZE_XPIC,xpicOption,NULL,timestep,0);
}

// vStar and vStarNext need to be zero on ghost nodes too
void XPICExtrapolationTask::InitializeXPICData(GridPatch *patchPtr,int xpicOption)
{	patchPtr->XPICSupport(INITIALIZE_XPIC,xpicOption,NULL,-1.,0);
}

// Double XPIC loop to find vStar
void XPICExtrapolationTask::XPICDoubleLoop(MPMBase *mpmptr,int matfld,int *nds,double *fn,int pn)
{
	// number of nodes
	int numnds = nds[0];

	// We building S^+Sv*(prev) on all nodes
	// The value for node i is Sum_p S_{ip}^+ S_{pj} v*(prev)_j
	// Here we are doing a single p and this is called for each particle
	// (Particle spin adds S_{ip}^{L+}} S_{pj}^{L} v*(prev)_j
	for(int i=1;i<=numnds;i++)
	{	// get mass from real node
		NodalPoint *ndptri = nd[nds[i]];
		short vfldi = mpmptr->vfld[i];
		double mass = ndptri->GetMaterialMass(vfldi,matfld);
		
		// get node pointer (which may now be ghost node)
		ndptri = GetNodePointer(pn,nds[i]);
		
		// loop over nodes again
		for(int j=1;j<=numnds;j++)
		{	// read-only real node to get vStarPrev
			NodalPoint *ndptrj = nd[nds[j]];
			short vfldj = mpmptr->vfld[j];
			Vector *vStarPrevj = ndptrj->GetVStarPrev(vfldj,matfld);

            // add to node i
            // Spi+ Spj = Mp Spi Spj/mi (Mp/mi applied at the end)
			double weight = fn[i]*fn[j];
			ndptri->AddVStarNext(vfldi,matfld,vStarPrevj,mpmptr->mp*weight/mass);;
		}
	}
}

// Transfer vStarNext to real node and zero it
void XPICExtrapolationTask::ReduceXPICData(GridPatch *patchPtr)
{	patchPtr->XPICSupport(COPY_VSTARNEXT,0,NULL,0.,0);
}

// Get Delta v on all nodes (vprev-nvext), set BC node velocities to zero,
//  and impose incrementl contact
void XPICExtrapolationTask::GetDeltaV(double timestep,int xpicOption)
{
	// Get Delta V =  (I-{S+}S)vprev = vprev-vnext (i.e., vnext = {S+}S vprev)
#pragma omp parallel for
	for(int i=1;i<=*nda;i++)
		nd[nda[i]]->XPICSupport(GET_DELTAV,0,NULL,timestep,0);
	
	// Set BC nodes to zero in latest Delta v
	NodalVelBC::GridVelocityConditions(xpicOption);
	
	// If contact (multimaterial or cracks) adjust contact calculations (mechanics only)
	ImposeIncrementalContact(timestep,xpicOption);
}

// Add Delta v (now in vprev to v(k) and zero vnext for next pass in the loop
// If contact nodes, add mass*Delta v to pk
void XPICExtrapolationTask::UpdateXStar(NodalPoint *ndptr,double timestep)
{	ndptr->XPICSupport(UPDATE_VSTAR,0,NULL,timestep,0);
}

// Change lastest increment in vk[VPREV] to be same in all contacting field
void XPICExtrapolationTask::ImposeIncrementalContact(double dtime,int callType)
{
	// recalculate interface energy each increment
	NodalPoint::interfaceEnergy=0.;
	
	// set incremental material contact nodes
	long numContactNodes = (long)MaterialContactNode::materialContactNodes.size();
	if(numContactNodes>0)
	{	int numNodesPerProc = (int)((double)numContactNodes/(double)(fmobj->GetNumberOfProcessors()));
#pragma omp parallel for if(numNodesPerProc>1)
		for(long i=0;i<numContactNodes;i++)
			MaterialContactNode::materialContactNodes[i]->NodalXPICIncrement(dtime,callType);
	}
	
	// set incremental crack contact nodes
	long numCrackNodes = (long)CrackNode::crackContactNodes.size();
	if(numCrackNodes>0)
	{	int numNodesPerProc = (int)((double)numCrackNodes/(double)(fmobj->GetNumberOfProcessors()));
#pragma omp parallel for if(numNodesPerProc>1)
		for(long i=0;i<numCrackNodes;i++)
			CrackNode::crackContactNodes[i]->NodalXPICIncrement(dtime,callType);
	}
}
