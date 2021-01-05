/********************************************************************************
	ExtrapolateRigidBCsTask.cpp
	nairn-mpm-fea

	Created by John Nairn on March 6, 2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.
 
	The tasks are:
	--------------
	* This task only active if there are rigid BC particles and the
	  simulation is to set them AFTER extrapolating to the grid
	  (these are from nmpmsRC to nmpms and excludes rigid block and contact particles)
	* First loop extrapolates rigid vel, T, and c to nodes
	* Second loop sets BCs on those nodes (and clears nodal values
	  because they are needed on subsequent steps and this step
	  is before any are needed).
	* May be velocity, temperature, and concentration
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/ExtrapolateRigidBCsTask.hpp"
#include "NairnMPM_Class/ProjectRigidBCsTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/RigidMaterial.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Boundary_Conditions/NodalTempBC.hpp"
#include "Boundary_Conditions/NodalConcBC.hpp"
#include "Custom_Tasks/ConductionTask.hpp"

#pragma mark CONSTRUCTORS

ExtrapolateRigidBCsTask::ExtrapolateRigidBCsTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get mass matrix, find dimensionless particle locations,
//	and find grid momenta
bool ExtrapolateRigidBCsTask::Execute(int taskOption)
{
#ifdef CONST_ARRAYS
	double fn[MAX_SHAPE_NODES];
#else
	double fn[maxShapeNodes];
#endif
	
	// undo dynamic velocity, temp, and conc BCs from rigid materials
	// and get pointer to first empty one in reuseRigid...BC
	ProjectRigidBCsTask::UnsetRigidBCs((BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
				  (BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
	ProjectRigidBCsTask::UnsetRigidBCs((BoundaryCondition **)&firstTempBC,(BoundaryCondition **)&lastTempBC,
				  (BoundaryCondition **)&firstRigidTempBC,(BoundaryCondition **)&reuseRigidTempBC);
	ProjectRigidBCsTask::UnsetRigidBCs((BoundaryCondition **)&firstConcBC,(BoundaryCondition **)&lastConcBC,
				  (BoundaryCondition **)&firstRigidConcBC,(BoundaryCondition **)&reuseRigidConcBC);
	
	int i,numnds,setFlags;
	Vector rvel;
	bool hasDir[3];
	double tempValue,concValue;
	
	// This loop not parallel because because it is normally small (only loops over rigid BC particles)
	// To make parallel, need ghost nodes and reduction step
	for(int p=nmpmsRC;p<nmpms;p++)
	{	// get material point and the rigid material
		MPMBase *mpmptr = mpm[p];
		const RigidMaterial *rigid = (RigidMaterial *)theMaterials[mpmptr->MatID()];				// material object for this particle
		
		// get rigid BC particle velocity
		// Get directions set, others will be zero
		// May not need to set to zero because setting BC will ignore those direction anyway
		setFlags = rigid->SetDirection();
		ZeroVector(&rvel);
		if(rigid->GetVectorSetting(&mpmptr->vel,hasDir,mtime,&mpmptr->pos))
		{	if(hasDir[0]) rvel.x = mpmptr->vel.x;
			if(hasDir[1]) rvel.y = mpmptr->vel.y;
			if(hasDir[2]) rvel.z = mpmptr->vel.z;
		}
		else
		{	if(setFlags&CONTROL_X_DIRECTION) rvel.x = mpmptr->vel.x;
			if(setFlags&CONTROL_Y_DIRECTION) rvel.y = mpmptr->vel.y;
			if(setFlags&CONTROL_Z_DIRECTION) rvel.z = mpmptr->vel.z;
		}
		
		// get rigid particle temperature
		if(rigid->RigidTemperature() && ConductionTask::active)
		{	setFlags += CONTROL_TEMPERATURE;
			if(rigid->GetValueSetting(&tempValue,mtime,&mpmptr->pos)) mpmptr->pTemperature = tempValue;
		}
		
		// concentration
		if(rigid->RigidConcentration() && fmobj->HasFluidTransport())
		{	setFlags += CONTROL_CONCENTRATION;
			if(rigid->GetValueSetting(&concValue,mtime,&mpmptr->pos)) mpmptr->pConcentration = concValue;
		}
		
		// get nodes and classic shape function for rigid material point p
		const ElementBase *elref = theElements[mpmptr->ElemID()];		// element containing this particle
		elref->ShapeFunction(mpmptr->GetNcpos(),false,&fn[1],NULL,NULL,NULL);
		numnds=elref->NumberNodes();
		
		// Add particle property to each node in the element
		for(i=1;i<=numnds;i++)
		{   // get node pointer and set values
			int mi=elref->nodes[i-1];		// 1 based node
			
			// it might possibly need a velocity field (does nothing in single material mode or if already there)
			nd[mi]->AddMatVelocityField(0,0);
			
			// add BC info
			nd[mi]->AddRigidBCInfo(mpmptr,fn[i],setFlags,&rvel);
		}
	}
	
	// read nodal settings, rezero their used values, and set boundary conditions
	for(int i=1;i<=nnodes;i++)
	{	NodalPoint *ndptr = nd[i];
		setFlags = ndptr->ReadAndZeroRigidBCInfo(&rvel,&tempValue,&concValue);
		if(setFlags==0) continue;

		if(setFlags&CONTROL_X_DIRECTION)
		{	ProjectRigidBCsTask::SetRigidBCs(i,-40,X_DIRECTION,rvel.x,0.,0,
						(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
						(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
		}
		if(setFlags&CONTROL_Y_DIRECTION)
		{	ProjectRigidBCsTask::SetRigidBCs(i,-40,Y_DIRECTION,rvel.y,0.,0,
						(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
						(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
		}
		if(setFlags&CONTROL_Z_DIRECTION)
		{	ProjectRigidBCsTask::SetRigidBCs(i,-40,Z_DIRECTION,rvel.z,0.,0,
						(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
						(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
		}
		if(setFlags&CONTROL_TEMPERATURE)
		{	ProjectRigidBCsTask::SetRigidBCs(i,0,TEMP_DIRECTION,tempValue,0.,0,
						(BoundaryCondition **)&firstTempBC,(BoundaryCondition **)&lastTempBC,
						(BoundaryCondition **)&firstRigidTempBC,(BoundaryCondition **)&reuseRigidTempBC);
		}
		if(setFlags&CONTROL_CONCENTRATION)
		{	ProjectRigidBCsTask::SetRigidBCs(i,0,CONC_DIRECTION,concValue,0.,0,
						(BoundaryCondition **)&firstConcBC,(BoundaryCondition **)&lastConcBC,
						(BoundaryCondition **)&firstRigidConcBC,(BoundaryCondition **)&reuseRigidConcBC);
		}
	}
	
	// if any left over rigid BCs, delete them now
	ProjectRigidBCsTask::RemoveRigidBCs((BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,(BoundaryCondition **)&firstRigidVelocityBC);
	ProjectRigidBCsTask::RemoveRigidBCs((BoundaryCondition **)&firstTempBC,(BoundaryCondition **)&lastTempBC,(BoundaryCondition **)&firstRigidTempBC);
	ProjectRigidBCsTask::RemoveRigidBCs((BoundaryCondition **)&firstConcBC,(BoundaryCondition **)&lastConcBC,(BoundaryCondition **)&firstRigidConcBC);

    return true;
}
