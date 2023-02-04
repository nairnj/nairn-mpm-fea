/********************************************************************************
	InitializationTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	The tasks are:
	--------------
	* Zero nodal values (InitialiseForTimeStep) for real and ghost nodes
		All MVF and CVF values
		Global transport values on node and for contact flow on CVF and MVF
	* Get shape function data (mostly CPDI, optionally for GIMP)
	* Set particle external forces
	* Clear out stored crack nodes and interface nodes
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/InitializationTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Exceptions/MPMWarnings.hpp"
#include "Boundary_Conditions/MatPtLoadBC.hpp"
#include "Cracks/CrackNode.hpp"
#include "Nodes/MaterialContactNode.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Patches/GridPatch.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark CONSTRUCTORS

InitializationTask::InitializationTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get mass matrix, find dimensionless particle locations,
//	and find grid momenta
// throws CommonException()
bool InitializationTask::Execute(int taskOption)
{
	CommonException *initErr = NULL;
	
	// Zero Mass Matrix and vectors
	warnings.BeginStep();
#pragma omp parallel
	{
		// zero active nodal variables on real nodes (first step does all)
		// After this step, nda is invalid until reset in post extrapolation task
#pragma omp for
		for(int i=1;i<=*nda;i++)
			nd[nda[i]]->InitializeForTimeStep();
		
        // zero ghost nodes in patch for this thread
        int pn = GetPatchNumber();
        patches[pn]->InitializeForTimeStep();
		
#pragma omp for nowait
		for(int p=0;p<nmpms;p++)
        {   MPMBase *mpmptr = mpm[p];                                       // pointer
			int elemIndex = mpmptr->ElemID();
			if(elemIndex==0) continue;
			const ElementBase *elref = theElements[elemIndex];		// element containing this particle
			try
			{	elref->GetShapeFunctionData(mpmptr);
			}
			catch(CommonException& err)
			{	if(initErr==NULL)
				{
#pragma omp critical (error)
					initErr = new CommonException(err);
				}
 			}
			catch(...)
			{	if(initErr==NULL)
				{
#pragma omp critical (error)
					initErr = new CommonException("Unexpected error","InitializationTask::Execute");
				}
 			}
		}
	}
	
	// was there an error?
	if(initErr!=NULL) throw *initErr;
    
    // Update forces applied to particles
	MatPtLoadBC::SetParticleFext(mtime);
	
	// remove nodes with contact
	if(fmobj->multiMaterialMode)
		MaterialContactNode::ReleaseContactNodes();
	if(firstCrack!=NULL)
		CrackNode::ReleaseContactNodes();
	
	// total interface energy
	NodalPoint::interfaceEnergy=0.;

    return true;
}
