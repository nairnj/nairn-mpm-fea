/********************************************************************************
	MassAndMomentumTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	This first MPM task projects mass and momentum to grid
 
	For each non-rigid particle:
		1. Update dilated volume on the particle (if needed)
		2. Extrapolate mass and momentum to the grid
		3. If needed extrapolate displacement and unscaled volume
		4. For multimaterial mode, extrapolate mass gradient
 
	For rigid multimaterial particles
		1. Find their current velocity
		2. Other the same (except mass treated different and transport is skipped)
 ********************************************************************************/

#include "NairnMPM_Class/MassAndMomentumTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/RigidMaterial.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Exceptions/CommonException.hpp"
#include "Patches/GridPatch.hpp"

#ifdef LOG_PROGRESS
#include "System/ArchiveData.hpp"
#endif

#pragma mark CONSTRUCTORS

MassAndMomentumTask::MassAndMomentumTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get mass matrix, find dimensionless particle locations,
//	and find grid momenta
// throws CommonException()
void MassAndMomentumTask::Execute(void)
{   
	CommonException *massErr = NULL;
    double fn[maxShapeNodes],xDeriv[maxShapeNodes],yDeriv[maxShapeNodes],zDeriv[maxShapeNodes];
    int ndsArray[maxShapeNodes];

	// loop over non-rigid and rigid contact particles - this parallel part changes only particle p
	// mass, momenta, etc are stored on ghost nodes, which are sent to real nodes in next non-parallel loop
    //for(int pn=0;pn<4;pn++)
#pragma omp parallel private(fn,xDeriv,yDeriv,zDeriv,ndsArray)
	{
        // thread for patch pn
		int pn = GetPatchNumber();
        
		// in case 2D planar
        for(int i=0;i<maxShapeNodes;i++) zDeriv[i] = 0.;
        
		try
		{	short vfld;
			NodalPoint *ndptr;
			int i,numnds,matfld,*nds;
			
			// Loop of non-rigid and rigid contact particles
			for(int block=FIRST_NONRIGID;block<=FIRST_RIGID_CONTACT;block++)
			{	MPMBase *mpmptr = patches[pn]->GetFirstBlockPointer(block);
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
						ndptr->AddMassMomentum(mpmptr,vfld,matfld,fn[i],xDeriv[i],yDeriv[i],zDeriv[i],
											   1,block==FIRST_NONRIGID);
					}
					
					// next material point
					mpmptr = (MPMBase *)mpmptr->GetNextObject();
				}
			}
		}
		catch(CommonException& err)
        {   if(massErr==NULL)
			{
#pragma omp critical (error)
				massErr = new CommonException(err);
			}
		}
		catch(...)
		{	if(massErr==NULL)
			{
#pragma omp critical (error)
				massErr = new CommonException("Unexpected error","MassAndMomentumTask::Execute");
			}
		}
	}
	
	// throw now - only possible error if too many CPDI nodes in 3D
	if(massErr!=NULL) throw *massErr;
    
	// reduction of ghost node forces to real nodes
	int totalPatches = fmobj->GetTotalNumberOfPatches();
	if(totalPatches>1)
	{	for(int pn=0;pn<totalPatches;pn++)
			patches[pn]->MassAndMomentumReduction();
	}
    
}

// Get Particle functions and constants
int MassAndMomentumTask::GetParticleFunctions(MPMBase *mpmptr,int **nds,double *fn,
											  double *xDeriv,double *yDeriv,double *zDeriv)
{
	const MaterialBase *matID = theMaterials[mpmptr->MatID()];		// material object for this particle
	int matfld = matID->GetField();										// material velocity field
		
	// get nodes and shape function for material point p
	const ElementBase *elref = theElements[mpmptr->ElemID()];
	if(fmobj->multiMaterialMode)
		elref->GetShapeGradients(fn,nds,xDeriv,yDeriv,zDeriv,mpmptr);
	else
		elref->GetShapeFunctions(fn,nds,mpmptr);
	
	return matfld;
}
