/********************************************************************************
	UpdateStrainsLastContactTask.cpp
	nairn-mpm-fea

	Created by John Nairn on 4/7/2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.

	Update strains on particles after they have been updated. This task
	is active for USL and USAVG update methods and when simulation
	requests a re=extrapolation to the grid (which often works better
	in problems with cracks or contact)
 
	The tasks are:
	--------------
	* For each non-rigid particle:
		- Extrapolate mass and momentum to the grid
		- If cracks or multimaterial mode
			+ Extrapolation position or displacement
			+ Extrapolate deformed area volume
			+ It multimaterial extrapolate volume gradient
	* Reduction to copy from ghost to real nodes
	* Update nodal transport values
	* Material contact and crack contact
	* Get grid BCs using GridMomentumConditions()
	* Full strain update
		- Get grid velocities (p/m)
		- Copy mechanical properties of a material
		- Tell material to update strain (etc) on the particle
		- Also extrapolates transport to particle and saves a "previous"
 
	This tasks not used when doing transport only
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/UpdateStrainsLastContactTask.hpp"
#include "NairnMPM_Class/UpdateStrainsFirstTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Custom_Tasks/TransportTask.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Cracks/CrackNode.hpp"
#include "Patches/GridPatch.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark CONSTRUCTORS

UpdateStrainsLastContactTask::UpdateStrainsLastContactTask(const char *name) : MassAndMomentumTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get total grid point forces (except external forces)
// throws CommonException()
void UpdateStrainsLastContactTask::Execute(void)
{
	CommonException *uslErr = NULL;
#ifdef CONST_ARRAYS
	int ndsArray[MAX_SHAPE_NODES];
	double fn[MAX_SHAPE_NODES],xDeriv[MAX_SHAPE_NODES],yDeriv[MAX_SHAPE_NODES],zDeriv[MAX_SHAPE_NODES];
#else
	int ndsArray[maxShapeNodes];
	double fn[maxShapeNodes],xDeriv[maxShapeNodes],yDeriv[maxShapeNodes],zDeriv[maxShapeNodes];
#endif
	
#pragma omp parallel private(ndsArray,fn,xDeriv,yDeriv,zDeriv)
	{
		// in case 2D planar
		for(int i=0;i<maxShapeNodes;i++) zDeriv[i] = 0.;
		
#pragma omp for
		// zero again (which finds new positions for contact rigid particle data on the nodes)
		for(int i=1;i<=nnodes;i++)
			nd[i]->RezeroNodeTask6(timestep);
		
		// zero ghost nodes on this patch
		int pn = GetPatchNumber();
		patches[pn]->RezeroNodeTask6(timestep);
		
		try
		{	short vfld;
			NodalPoint *ndptr;
			int i,numnds,matfld,*nds;
			
			// loop over non-rigid particles only
			MPMBase *mpmptr = patches[pn]->GetFirstBlockPointer(FIRST_NONRIGID);
			while(mpmptr!=NULL)
			{	// get shape functions and mat field
				nds = ndsArray;
				matfld = GetParticleFunctions(mpmptr,&nds,fn,xDeriv,yDeriv,zDeriv);
				numnds = nds[0];
				
				// Add particle property to each node in the element
				for(i=1;i<=numnds;i++)
				{   // get node pointer
					ndptr = GetNodePointer(pn,nds[i]);
					
					// add mass and momentum (and maybe contact stuff) to this node
					vfld = mpmptr->vfld[i];
					ndptr->AddMassMomentumLast(mpmptr,vfld,matfld,fn[i],xDeriv[i],yDeriv[i],zDeriv[i]);
				}
				
				// next material point
				mpmptr = (MPMBase *)mpmptr->GetNextObject();
			}
		}
		catch(CommonException& err)
		{	if(uslErr==NULL)
			{
#pragma omp critical (error)
				uslErr = new CommonException(err);
			}
		}
		catch(std::bad_alloc&)
		{	if(uslErr==NULL)
			{
#pragma omp critical (error)
				uslErr = new CommonException("Memory error","UpdateStrainsLastContactTask::Execute");
			}
		}
		catch(...)
		{	if(uslErr==NULL)
			{
#pragma omp critical (error)
				uslErr = new CommonException("Unexpected error","UpdateStrainsLastContactTask::Execute");
			}
		}
	}
	
	// throw errors now
	if(uslErr!=NULL) throw *uslErr;
	
	// reduction of ghost node forces to real nodes
	int totalPatches = fmobj->GetTotalNumberOfPatches();
	if(totalPatches>1)
	{	for(int pn=0;pn<totalPatches;pn++)
			patches[pn]->MassAndMomentumReductionLast();
	}
	
	// grid temperature is never updated unless needed here
	// update nodal values for transport properties (when coupled to strain)
	TransportTask *nextTransport=transportTasks;
	while(nextTransport!=NULL)
		nextTransport=nextTransport->UpdateNodalValues(timestep);
	
	// adjust momenta for multimaterial contact
	if(fmobj->multiMaterialMode)
	{	for(int i=1;i<=nnodes;i++)
			nd[i]->MaterialContactOnNode(timestep,UPDATE_STRAINS_LAST_CALL);
	}
	
	// adjust momenta for crack contact
	if(firstCrack!=NULL) CrackNode::ContactOnKnownNodes();
	
	// impose grid boundary conditions
	NodalVelBC::GridMomentumConditions();
	
	// update strains based on current velocities
	UpdateStrainsFirstTask::FullStrainUpdate(strainTimestep,(fmobj->mpmApproach==USAVG_METHOD),fmobj->np);
}
