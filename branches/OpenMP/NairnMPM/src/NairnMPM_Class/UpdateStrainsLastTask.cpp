/********************************************************************************
	UpdateStrainsLastTask.cpp
	NairnMPM

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	A strain update at the end of the MPM step is used in the SZS and
		the USAVG method. For these to be stable, the new particle momenta
		have to be re-extrapolated to the grid
 
	Zero nodes and extrpolate particle momenta and and positions
		to the grid.
 
	Once done, adject momenta for material contact, crack contact, and
		boundary conditions
	
	Since will reuse initial locations (i.e. shape functions and their gradients
		on each particle are the same), the mass, volume, and mass gradient
		would not change, even if re-extrapolated, and thus they do not
		need to be changed.
 
	For rigid particles, only displacement changed and it is found when node
		rezeroed
 
	After new extrapolations, update strains on all particles
********************************************************************************/

#include "NairnMPM_Class/UpdateStrainsLastTask.hpp"
#include "NairnMPM_Class/UpdateStrainsFirstTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Custom_Tasks/TransportTask.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Cracks/CrackNode.hpp"

// NEWINCLUDE
#include "Patches/GridPatch.hpp"

#pragma mark CONSTRUCTORS

UpdateStrainsLastTask::UpdateStrainsLastTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get total grid point forces (except external forces)
void UpdateStrainsLastTask::Execute(void)
{
	int nds[maxShapeNodes];
	double fn[maxShapeNodes],xDeriv[maxShapeNodes],yDeriv[maxShapeNodes],zDeriv[maxShapeNodes];
	
#pragma omp parallel private(nds,fn,xDeriv,yDeriv,zDeriv)
	{
#pragma omp for
        // zero again (which finds new positions for rigid particles)
		for(int i=1;i<=nnodes;i++)
			nd[i]->RezeroNodeTask6(timestep);
		
#ifdef _OPENMP
		int pn = omp_get_thread_num();
#else
		int pn = 0;
#endif
/*
	int tp = fmobj->GetTotalNumberOfPatches();
	for(int pn=0;pn<tp;pn++)
	{
*/
        // zero ghost nodes on this patch
		patches[pn]->RezeroNodeTask6(timestep);
        
        // loop over non-rigid particles only - this parallel part changes only particle p
        // mass, momenta, etc are stored on ghost nodes, which are sent to real nodes in next non-parallel loop
        MPMBase *mpmptr = patches[pn]->GetFirstBlockPointer(FIRST_NONRIGID);
        while(mpmptr!=NULL)
        {   const MaterialBase *matref = theMaterials[mpmptr->MatID()];
            int matfld = matref->GetField();
            
            // find shape functions (why ever need gradients?)
            const ElementBase *elref = theElements[mpmptr->ElemID()];
            int numnds;
            if(fmobj->multiMaterialMode)
            {   // Need gradients for volume gradient
                elref->GetShapeGradients(&numnds,fn,nds,mpmptr->GetNcpos(),xDeriv,yDeriv,zDeriv,mpmptr);
            }
            else
                elref->GetShapeFunctions(&numnds,fn,nds,mpmptr->GetNcpos(),mpmptr);
            
            short vfld;
            NodalPoint *ndptr;
            for(int i=1;i<=numnds;i++)
            {
#ifdef _OPENMP
                ndptr = patches[pn]->GetNodePointer(nds[i]);
#else
                ndptr = nd[nds[i]];
#endif
                vfld = (short)mpmptr->vfld[i];
                ndptr->AddMassMomentumLast(mpmptr,vfld,matfld,fn[i],xDeriv[i],yDeriv[i],zDeriv[i]);

            }
            
            // next non-rigid material point
            mpmptr = (MPMBase *)mpmptr->GetNextObject();
        }
	}
	
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
			nd[i]->MaterialContactOnNode(timestep,TRUE,NULL,NULL);
	}
	
	// adjust momenta for crack contact
	if(firstCrack!=NULL) CrackNode::ContactOnKnownNodes();
	
	// impose grid boundary conditions
	NodalVelBC::GridMomentumConditions(FALSE);
	
	// update strains based on current velocities
	UpdateStrainsFirstTask::FullStrainUpdate(strainTimestep,(fmobj->mpmApproach==USAVG_METHOD),fmobj->np);
}
