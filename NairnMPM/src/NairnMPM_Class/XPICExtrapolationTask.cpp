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

// class globals
XPICExtrapolationTask *XPICMechanicsTask=NULL;

#pragma mark CONSTRUCTORS

XPICExtrapolationTask::XPICExtrapolationTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Update particle position, velocity, temp, and conc
// throws CommonException()
bool XPICExtrapolationTask::Execute(int xpicOption)
{
	// get order for this version of XPIC and skip if no calculations needed
	int m = GetXPICOrder();
	if(m==0) return true;
	
	CommonException *xpicErr = NULL;
#ifdef CONST_ARRAYS
	double fn[MAX_SHAPE_NODES],xDeriv[MAX_SHAPE_NODES],yDeriv[MAX_SHAPE_NODES],zDeriv[MAX_SHAPE_NODES];
	int ndsArray[MAX_SHAPE_NODES];
#else
	double fn[maxShapeNodes],xDeriv[maxShapeNodes],yDeriv[maxShapeNodes],zDeriv[maxShapeNodes];
	int ndsArray[maxShapeNodes];
#endif
	
	// zero vStar and vStarNext on real nodes and
	// initialize vStarPrev on MVF for each real node to pi0/mi = (pi-fi*dt)/mi
#pragma omp parallel for
	for(int i=1;i<=*nda;i++)
		InitializeXPICData(nd[nda[i]],timestep,xpicOption);
	
	// In parallel, initialize patches too
	int totalPatches = fmobj->GetTotalNumberOfPatches();
	if(totalPatches>1)
	{	for(int i=0;i<totalPatches;i++)
			InitializeXPICData(patches[i],xpicOption);
	}
	
	// iterate for k from 2 to XPIC order m (note that loop is skipped for order 1)
	double vsign = -1.;					// (-1)^k starting at -1 for k=2 to subtract v*
	for(int k=2;k<=m;k++)
	{
#pragma omp parallel private(fn,ndsArray,xDeriv,yDeriv,zDeriv)
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
					if(XPICDoubleLoopNeedsGradients())
						elref->GetShapeGradients(fn,&nds,xDeriv,yDeriv,zDeriv,mpmptr);
					else
						elref->GetShapeFunctions(fn,&nds,mpmptr);
					
					// Add to each node for this particle
					double scale = (double)(m-k+1)/(double)k;
#if MM_XPIC == 1
					double scaleContact = (double)(m-k)/(double)k;
#else
					double scaleContact = 1.;
#endif
					
					// double loop over nodes
					XPICDoubleLoop(mpmptr,matfld,nds,fn,pn,scale,scaleContact,xDeriv,yDeriv,zDeriv);
					
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
		
		// reduction of ghost node velocities to real nodes (and zero vStarNext on ghost nodes)
		if(totalPatches>1)
		{	for(int pn=0;pn<totalPatches;pn++)
				ReduceXPICData(patches[pn],k);
		}
		
		// Increment vStar and copy vStarNext to vStarPrev (which is only needed on real nodes) and zero vStarNext
#pragma omp parallel for
		for(int i=1;i<=*nda;i++)
			UpdateXStar(nd[nda[i]],timestep,m,k,vsign);
		
		// change the sign
		vsign = -vsign;
	}
	
	// Done unless this XPIC wants to extrapolated back to particles
	if(!XPICDoesBackExtrapolation()) return true;
	
	// Extrapolate back tothe particles
#pragma omp parallel private(fn,ndsArray)
	{
		// thread for patch pn
		int pn = GetPatchNumber();
		
		try
		{	// Loop over non-rigid particles
			MPMBase *mpmptr = patches[pn]->GetFirstBlockPointer(FIRST_NONRIGID);
			while(mpmptr!=NULL)
			{	// get nodes and shape function for material point p
				const ElementBase *elref = theElements[mpmptr->ElemID()];
				int *nds = ndsArray;
				elref->GetShapeFunctions(fn,&nds,mpmptr);
				
				// Back Extrapolations
				XPICBackExtrapolation(mpmptr,nds,fn,m);
				
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
	
	// throw now - only possible error if too many CPDI nodes in 3D
	if(xpicErr!=NULL) throw *xpicErr;

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
{	ndptr->XPICSupport(INITIALIZE_XPIC,xpicOption,NULL,timestep,0,0,0.);
}

// vStar and vStarNext need to be zero on ghost nodes too
void XPICExtrapolationTask::InitializeXPICData(GridPatch *patchPtr,int xpicOption)
{	patchPtr->XPICSupport(INITIALIZE_XPIC,xpicOption,NULL,-1.,0,0,0.);
}

// need gradients in XPIC double loop
bool XPICExtrapolationTask::XPICDoubleLoopNeedsGradients(void) { return true; }

// Double XPIC loop to find vStar
void XPICExtrapolationTask::XPICDoubleLoop(MPMBase *mpmptr,int matfld,int *nds,double *fn,int pn,double scale,
										   double scaleContact,double *xDeriv,double *yDeriv,double *zDeriv)
{
	// number of nodes
	int numnds = nds[0];
	
	for(int i=1;i<=numnds;i++)
	{	// get mass from real node
		NodalPoint *ndptri = nd[nds[i]];
		short vfldi = mpmptr->vfld[i];
		double mass = ndptri->GetMaterialMass(vfldi,matfld);
		Vector *delXiMpPtr = NULL;
		
		// get node pointer (which may now be ghost node)
		ndptri = GetNodePointer(pn,nds[i]);
		
		// loop over nodes again
		for(int j=1;j<=numnds;j++)
		{	// read-only real node to get vStarPrev
			NodalPoint *ndptrj = nd[nds[j]];
			short vfldj = mpmptr->vfld[j];
			Vector *vStarPrevj = ndptrj->GetVStarPrev(vfldj,matfld);
			Vector *delXjPtr = NULL;

            // add to node i
			double weight = mpmptr->mp*fn[i]*fn[j]/mass;
			ndptri->AddVStarNext(vfldi,matfld,vStarPrevj,delXiMpPtr,delXjPtr,NULL,scale*weight,scaleContact*weight);
		}
	}
}

// Transfer vStarNext to real node and zero it
void XPICExtrapolationTask::ReduceXPICData(GridPatch *patchPtr,int k)
{	patchPtr->XPICSupport(COPY_VSTARNEXT,0,NULL,0.,0,k,0.);
}

// Update vStar, transfer vStarNext to vStarPrev and zero vStarNext
void XPICExtrapolationTask::UpdateXStar(NodalPoint *ndptr,double timestep,int m,int k,double vsign)
{	ndptr->XPICSupport(UPDATE_VSTAR,0,NULL,timestep,m,k,vsign);
}

// Does not due back extrapolate
bool XPICExtrapolationTask::XPICDoesBackExtrapolation(void)
{	return false;
}

// empty method
void XPICExtrapolationTask::XPICBackExtrapolation(MPMBase *mpmptr,int *nds,double *fn,int m) {}

