/********************************************************************************
	InitializationTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
********************************************************************************/

#include "NairnMPM_Class/InitializationTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Nodes/MaterialInterfaceNode.hpp"
#include "Exceptions/MPMWarnings.hpp"
#include "Boundary_Conditions/MatPtLoadBC.hpp"
#include "Cracks/CrackNode.hpp"
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
void InitializationTask::Execute(void)
{
	CommonException *initErr = NULL;
	
	// Zero Mass Matrix and vectors
	warnings.BeginStep();
	
#pragma omp parallel
	{
        // zero all nodal variables on real nodes
#pragma omp for
		for(int i=1;i<=nnodes;i++)
			nd[i]->InitializeForTimeStep();
		
        // zero ghost nodes in patch for this thread
        int pn = GetPatchNumber();
        patches[pn]->InitializeForTimeStep();

		// particle calculations get CPDI or GIMP info for each nonrigid and rigid contact particle
#pragma omp for nowait
		for(int p=0;p<nmpmsRC;p++)
        {   MPMBase *mpmptr = mpm[p];                                       // pointer
			const ElementBase *elref = theElements[mpmptr->ElemID()];		// element containing this particle
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
	
	// remove contact conditions
	CrackNode::RemoveCrackNodes();
	MaterialInterfaceNode::RemoveInterfaceNodes();
}	
