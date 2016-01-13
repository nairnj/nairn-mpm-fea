/******************************************************************************************
	GridForcesTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	Find all forces on the grid including internal forces (from particle stress)
		external forces, and body forces. Add equivalent forces for transport
		tasks.
 
	After main loop, get traction forces on cracks, add crack tip heating,
		and imperfect interface forces on cracks and track interface energy.
 
	Sum all forces only with external damping on the nodes
 
	Finally reconcile forces with boundary conditions. Do same for transport
		tasks.
******************************************************************************************/

#include "NairnMPM_Class/GridForcesTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/TransportTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Patches/GridPatch.hpp"
#include "Exceptions/CommonException.hpp"
#ifdef LOG_PROGRESS
#include "System/ArchiveData.hpp"
#endif

#pragma mark CONSTRUCTORS

GridForcesTask::GridForcesTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get total grid point forces (except external forces)
void GridForcesTask::Execute(void)
{
	CommonException *forceErr = NULL;
	
	// need to be private in threads
	TransportProperties t;
	int numnds,nds[maxShapeNodes];
	double fn[maxShapeNodes],xDeriv[maxShapeNodes],yDeriv[maxShapeNodes],zDeriv[maxShapeNodes];
	
	// loop over non-rigid particles - this parallel part changes only particle p
	// forces are stored on ghost nodes, which are sent to real nodes in next non-parallel loop
#pragma omp parallel private(t,numnds,nds,fn,xDeriv,yDeriv,zDeriv)
	{	// in case 2D planar
        for(int i=0;i<maxShapeNodes;i++) zDeriv[i] = 0.;
        
        // patch for this thread
        int pn = GetPatchNumber();
        
		try
		{	MPMBase *mpmptr = patches[pn]->GetFirstBlockPointer(FIRST_NONRIGID);
			while(mpmptr!=NULL)
			{	const MaterialBase *matref = theMaterials[mpmptr->MatID()];		// material class (read only)
				int matfld = matref->GetField(); 
				
				// get transport tensors (if needed)
				if(transportTasks!=NULL)
					matref->GetTransportProps(mpmptr,fmobj->np,&t);
				
				// find shape functions and derviatives
				const ElementBase *elemref = theElements[mpmptr->ElemID()];
				elemref->GetShapeGradients(&numnds,fn,nds,xDeriv,yDeriv,zDeriv,mpmptr);
				
				// Add particle property to buffer on the material point (needed to allow parallel code)
				short vfld;
				NodalPoint *ndptr;
				for(int i=1;i<=numnds;i++)
				{	vfld = (short)mpmptr->vfld[i];					// crack velocity field to use
					
					// total force vector = internal + external forces
					//	(in g mm/sec^2 or micro N)
					Vector theFrc;
					mpmptr->GetFintPlusFext(&theFrc,fn[i],xDeriv[i],yDeriv[i],zDeriv[i]);
					
					// add body forces (do in outside loop now)
					
					// add the total force to nodal point
                    ndptr = GetNodePointer(pn,nds[i]);
					ndptr->AddFtotTask3(vfld,matfld,&theFrc);
					
#ifdef CHECK_NAN
                    if(theFrc.x!=theFrc.x || theFrc.y!=theFrc.y || theFrc.z!=theFrc.z)
                    {
#pragma omp critical (output)
						{	cout << "\n# GridForcesTask::Execute: bad nodal force vfld = " << vfld << ", matfld = " << matfld;
							PrintVector(" theFrc = ",&theFrc);
							cout << endl;
							ndptr->Describe();
						}
                    }
#endif
					// transport forces
					TransportTask *nextTransport=transportTasks;
					while(nextTransport!=NULL)
						nextTransport=nextTransport->AddForces(ndptr,mpmptr,fn[i],xDeriv[i],yDeriv[i],zDeriv[i],&t);
				}

				// next material point
				mpmptr = (MPMBase *)mpmptr->GetNextObject();
			}
		}
		catch(CommonException err)
		{	if(forceErr==NULL)
			{
#pragma omp critical (error)
				forceErr = new CommonException(err);
			}
		}
	}
	
	// throw errors now
	if(forceErr!=NULL) throw *forceErr;
	
	// reduction of ghost node forces to real nodes
	int totalPatches = fmobj->GetTotalNumberOfPatches();
	if(totalPatches>1)
	{	for(int pn=0;pn<totalPatches;pn++)
			patches[pn]->GridForcesReduction();
	}
	
}
