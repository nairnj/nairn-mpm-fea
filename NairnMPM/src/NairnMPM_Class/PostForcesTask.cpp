/********************************************************************************
	PostForcesTask.cpp
	nairn-mpm-fea

	Created by John Nairn on March 8, 2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.

	Tasks done after extrapolating forces to the grid:
	-------------------------------------------------
	* Add traction BC nodal forces
	* Add crack tractions to nodal forces
	* Add crack tip heating to conduction force
	* Add imperfectinterface forces to nodes
	* Add gravity and body forces
	* Make nodes with velocity BCs have consistent forces
	* Make transport tasks forces consistent with transport nodal BCs
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/PostForcesTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "Boundary_Conditions/MatPtTractionBC.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Exceptions/CommonException.hpp"
#ifdef RESTART_OPTION
#include "NairnMPM_Class/MeshInfo.hpp"
#endif

#pragma mark CONSTRUCTORS

PostForcesTask::PostForcesTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get mass matrix, find dimensionless particle locations,
//	and find grid momenta
bool PostForcesTask::Execute(int taskOption)
{
	// restore nodal momenta
#pragma omp parallel for
	for(int i=1;i<=*nda;i++)
		nd[nda[i]]->RestoreMomenta();
	
	// Add traction BCs on particles
	MatPtTractionBC::SetParticleSurfaceTractions(mtime);
	
	// Add traction law forces to velocity fields
	if(fmobj->hasTractionCracks)
	{	CrackHeader *nextCrack=firstCrack;
		while(nextCrack!=NULL)
		{	nextCrack->AddTractionForce();
			nextCrack=(CrackHeader *)nextCrack->GetNextObject();
		}
	}
	
	// Add crack tip heating adds to conduction force
	if(conduction) conduction->AddCrackTipHeating();

	// Add gravity and body forces (if any are present)
	// Note: If ever need to implement body force that depend on particle state (stress, strain, etc.)
	//			then move the body force addition into GridForcesTask loop where gravity is commented out
	Vector gridBodyForce;
	if(bodyFrc.HasGridDampingForces())
	{	CommonException *bfErr = NULL;

#pragma omp parallel for private(gridBodyForce)
		for(int i=1;i<=*nda;i++)
		{	NodalPoint *ndptr = nd[nda[i]];
			try
			{	Vector fpos = MakeVector(ndptr->x,ndptr->y,ndptr->z);
				bodyFrc.GetGridBodyForce(&gridBodyForce,&fpos,mtime);
				ndptr->AddGravityAndBodyForceTask3(&gridBodyForce);
			}
			catch(CommonException &err)
			{   if(bfErr==NULL)
				{
#pragma omp critical (error)
					bfErr = new CommonException(err);
				}
 			}
		}
			
		// throw now - only known error is problem with function for body force setting
		if(bfErr!=NULL) throw *bfErr;
	}

    // Impose BCs on ftot to get correct grid BCs for velocity
	if(bodyFrc.GetXPICOrder()<=1 || bodyFrc.GridBCOption()!=GRIDBC_VELOCITY_ONLY)
		NodalVelBC::GridVelocityConditions(GRID_FORCES_CALL);

	// Transport force BCs
	TransportTask::TransportForceBCs(timestep);

#ifdef RESTART_OPTION
    if(fabs(fmobj->restartScaling)>1.e-6)
    {   // check if any node is moving too fast
        double maxDist = fmobj->restartCFL*mpmgrid.GetMinCellDimension();
        for(int i=1;i<=*nda;i++)
        {   if(nd[nda[i]]->IsTravelTooMuch(timestep,maxDist))
                return false;
        }
    }
#endif
    
    return true;
}
