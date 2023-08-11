/********************************************************************************
	ProjectRigidBCsTask.cpp
	nairn-mpm-fea
 
	Created by John Nairn on March 6, 2015
	Copyright (c) 2015 John A. Nairn, All rights reserved.
 
	The tasks are:
	-------------
	* This task only active when have some rigid BC particles
			(particles from nmpmsRC to nmpms-1)
	* Project values of rigid particles to all nodes in the element
	  containing the particle
	* May be velocity, temperature, and concentration
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/ProjectRigidBCsTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/RigidMaterial.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Boundary_Conditions/NodalTempBC.hpp"
#include "Boundary_Conditions/NodalConcBC.hpp"
#include "Exceptions/CommonException.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"

#pragma mark CONSTRUCTORS

ProjectRigidBCsTask::ProjectRigidBCsTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get mass matrix, find dimensionless particle locations,
//	and find grid momenta
bool ProjectRigidBCsTask::Execute(int taskOption)
{
#ifdef CONST_ARRAYS
	double fn[MAX_SHAPE_NODES];
	int ndsArray[MAX_SHAPE_NODES];
#else
	double fn[maxShapeNodes];
	int ndsArray[maxShapeNodes];
#endif
	
	// undo dynamic velocity, temp, and conc BCs from rigid materials
	// and get pointer to first empty one in reuseRigid...BC
	UnsetRigidBCs((BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
				  (BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
	UnsetRigidBCs((BoundaryCondition **)&firstTempBC,(BoundaryCondition **)&lastTempBC,
				  (BoundaryCondition **)&firstRigidTempBC,(BoundaryCondition **)&reuseRigidTempBC);
	UnsetRigidBCs((BoundaryCondition **)&firstConcBC,(BoundaryCondition **)&lastConcBC,
				  (BoundaryCondition **)&firstRigidConcBC,(BoundaryCondition **)&reuseRigidConcBC);

	// For Rigid BC materials create velocity BC on each node in the element
	for(int p=nmpmsRC;p<nmpms;p++)
	{	MPMBase *mpmptr = mpm[p];										// pointer
		if(mpmptr->InReservoir()) continue;
		
		int matid0 = mpmptr->MatID();
		const MaterialBase *matID = theMaterials[matid0];				// material object for this particle
		RigidMaterial *rigid=(RigidMaterial *)matID;

		const ElementBase *elref = theElements[mpmptr->ElemID()];		// element containing this particle
		int *nds = ndsArray;
		elref->GetShapeFunctions(fn,&nds,mpmptr);
		int numnds = nds[0];

		double rvalue;
		for(int i=1;i<=numnds;i++)
		{
			int mi=nds[i];
			
			// look for setting function in one to three directions in rigid BC particle
			// GetVectorSetting() returns true if function has set the velocity, otherwise it returns FALSE
			bool hasDir[3];
			if(rigid->GetVectorSetting(&mpmptr->vel,hasDir,mtime,&mpmptr->pos))
			{   // velocity set by 1 to 3 functions as determined by hasDir[i]
				if(hasDir[0])
				{	SetRigidBCs(mi,matid0,X_DIRECTION,mpmptr->vel.x,0.,rigid->mirrored,
								(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
								(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
				}
				if(hasDir[1])
				{	SetRigidBCs(mi,matid0,Y_DIRECTION,mpmptr->vel.y,0.,rigid->mirrored,
								(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
								(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
				}
				if(hasDir[2])
				{	SetRigidBCs(mi,matid0,Z_DIRECTION,mpmptr->vel.z,0.,rigid->mirrored,
								(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
								(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
				}
			}
			else
			{   // velocity set by particle velocity in selected directions
				if(rigid->RigidDirection(X_DIRECTION))
				{	SetRigidBCs(mi,matid0,X_DIRECTION,mpmptr->vel.x,0.,rigid->mirrored,
								(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
								(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
				}
				if(rigid->RigidDirection(Y_DIRECTION))
				{	SetRigidBCs(mi,matid0,Y_DIRECTION,mpmptr->vel.y,0.,rigid->mirrored,
								(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
								(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
				}
				if(rigid->RigidDirection(Z_DIRECTION))
				{	SetRigidBCs(mi,matid0,Z_DIRECTION,mpmptr->vel.z,0.,rigid->mirrored,
								(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
								(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
				}
			}
			
			// temperature
			if(rigid->RigidTemperature())
			{   // if has function change it, otherwise leave as is
                if(rigid->GetValueSetting(&rvalue,mtime,&mpmptr->pos)) mpmptr->pTemperature=rvalue;
				SetRigidBCs(mi,matid0,TEMP_DIRECTION,mpmptr->pTemperature,0.,0,
							(BoundaryCondition **)&firstTempBC,(BoundaryCondition **)&lastTempBC,
							(BoundaryCondition **)&firstRigidTempBC,(BoundaryCondition **)&reuseRigidTempBC);
			}
			
			// concentration (but only if has diffusion or poroelasticity)
            // only concentratation and always in pDiff[0]
			if(rigid->RigidConcentration())
			{   // if has function change it, otherwise leave as is
                if(rigid->GetValueSetting(&rvalue,mtime,&mpmptr->pos)) mpmptr->pDiff[0]->conc=rvalue;
				if(fmobj->HasDiffusion())
				{	if(mpmptr->pDiff[0]->conc<0.)
						mpmptr->pDiff[0]->conc=0.;
					else if(!diffusion->noLimit && mpmptr->pDiff[0]->conc>1.)
						mpmptr->pDiff[0]->conc=1.;
				}
 				SetRigidBCs(mi,matid0,CONC_DIRECTION,mpmptr->pDiff[0]->conc,0.,0,
							(BoundaryCondition **)&firstConcBC,(BoundaryCondition **)&lastConcBC,
							(BoundaryCondition **)&firstRigidConcBC,(BoundaryCondition **)&reuseRigidConcBC);
			}
		}
	}
	
	// if any left over rigid BCs, delete them now
	RemoveRigidBCs((BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,(BoundaryCondition **)&firstRigidVelocityBC);
	RemoveRigidBCs((BoundaryCondition **)&firstTempBC,(BoundaryCondition **)&lastTempBC,(BoundaryCondition **)&firstRigidTempBC);
	RemoveRigidBCs((BoundaryCondition **)&firstConcBC,(BoundaryCondition **)&lastConcBC,(BoundaryCondition **)&firstRigidConcBC);
    
    // reset pointer to first one
    if(fmobj->HasFluidTransport())
    {   firstDiffBC[MOISTURE_DIFFUSION] = firstConcBC;
    }
    else
    {   firstDiffBC[POROELASTICITY_DIFFUSION] = firstConcBC;
    }

    return true;
}

// Unset nodal dof for rigid BCs so can try to reuse
//   them without needing new memory allocations
// No rigid BCs are deleted and when done, reuseRigidBC is set to first empty one (excess deleted later)
// Note that set BC will set a direction thereby preventing any dynamic BC in the same direction
//		Thus dynamic ones should not be in the same direction
//		Furthermore, nodal set direction will always be correct indicator of any BC (set of dynamic) on that node
void ProjectRigidBCsTask::UnsetRigidBCs(BoundaryCondition **firstBC,BoundaryCondition **lastBC,
										BoundaryCondition **firstRigidBC,BoundaryCondition **reuseRigidBC)
{
	// exit if none
	if(*firstRigidBC==NULL) return;
	
	// unset dynamic ones
    BoundaryCondition *nextBC=*firstRigidBC;
	while(nextBC!=NULL)
		nextBC=nextBC->UnsetDirection();
	
	// were they all rigid BCs, but still has some?
	if(*firstBC==*firstRigidBC)
	{	*lastBC=NULL;				// since none set yet
		*reuseRigidBC=*firstRigidBC;
		return;
	}
	
	// otherwise search for last actual grid BC
    nextBC=*firstBC;
	while(nextBC->GetNextObject()!=*firstRigidBC)
		nextBC=(BoundaryCondition *)nextBC->GetNextObject();
	*lastBC=nextBC;
	*reuseRigidBC=*firstRigidBC;
}

// Set boundary conditions determined by moving rigid paticles
// throws std::bad_alloc
void ProjectRigidBCsTask::SetRigidBCs(int mi,int matid0,int type,double value,double angle,int mirrored,
									  BoundaryCondition **firstBC,BoundaryCondition **lastBC,
									  BoundaryCondition **firstRigidBC,BoundaryCondition **reuseRigidBC)
{
	BoundaryCondition *newBC=NULL;
	
	// check if already set in that direction by actual BC or by previous rigid BC
	// New rigid BC's can only be on free directions
	if(nd[mi]->fixedDirection&type) return;
	
	// create new boundary conditions
	switch(type)
    {   case X_DIRECTION:
		case Y_DIRECTION:
		case Z_DIRECTION:
			if(*reuseRigidBC!=NULL)
				newBC=(*reuseRigidBC)->SetRigidProperties(mi,type,CONSTANT_VALUE,value);
			else
            {   int newType = type==Z_DIRECTION ? Z_DIRECTION_INPUT : type ;
				newBC=(BoundaryCondition *)(new NodalVelBC(mi,newType,CONSTANT_VALUE,value,(double)0.,(double)0.,(double)0.));
			}
			((NodalVelBC *)newBC)->SetMirrorSpacing(mirrored);
			break;
			
		case TEMP_DIRECTION:
			if(*reuseRigidBC!=NULL)
				newBC=(*reuseRigidBC)->SetRigidProperties(mi,type,CONSTANT_VALUE,value);
			else
            {   newBC=(BoundaryCondition *)(new NodalTempBC(mi,CONSTANT_VALUE,value,(double)0.));
            }
			break;
			
		case CONC_DIRECTION:
 			if(*reuseRigidBC!=NULL)
				newBC=(*reuseRigidBC)->SetRigidProperties(mi,type,CONSTANT_VALUE,value);
			else
            {   newBC=(BoundaryCondition *)(new NodalConcBC(mi,CONSTANT_VALUE,value,
                                                            (double)0.,MOISTURE_DIFFUSION));
            }
			break;
			
		default:
			break;
			
	}
	
	// set BC id in case needed to 1-based material number
	newBC->SetID(matid0+1);
	
	// *firstBC and *lastBC are first and last of this type
	// *firstRigidBC will save the first one
	// if *reuseRigidBC!=NULL, then reusing previous rigid BCs
	if(*firstBC==NULL)
	{	// Only happens when no normal BCs and no rigidBCs to reuse so start with this one
		*firstBC=newBC;
		*firstRigidBC=newBC;
		// reuseRigidBC must be NULL already
	}
	else
	{	if(*reuseRigidBC!=NULL)
		{	// next object of last BC is already set
			// firstRigidBC is already valid
			// advance to reuse next one (or could get to NULL if all used up)
			*reuseRigidBC=(BoundaryCondition *)(*reuseRigidBC)->GetNextObject();
		}
		else
		{	// created a new one or ran out of ones to reuse
			(*lastBC)->SetNextObject(newBC);
			if(*firstRigidBC==NULL) *firstRigidBC=newBC;
		}
	}
	*lastBC=newBC;
}

// Remove any dynamically created boundary conditions
//    that are no longer needed
void ProjectRigidBCsTask::RemoveRigidBCs(BoundaryCondition **firstBC,BoundaryCondition **lastBC,BoundaryCondition **firstRigidBC)
{
	// exit if none
	if(*firstRigidBC==NULL) return;
	
	// trap if this step did not reuse any BCs
	BoundaryCondition *nextBC,*prevBC;
	if(*lastBC==NULL)
	{	nextBC=*firstBC;		// will be deleting them all below
		// None reused and no grid ones either?
		*firstBC=NULL;
		*firstRigidBC=NULL;
	}
	else
	{	nextBC=(BoundaryCondition *)(*lastBC)->GetNextObject();
		(*lastBC)->SetNextObject(NULL);			// next one on lastBC needs to be NULL
		// Were none of the rigid ones resued? (they are deleted below)
		if(nextBC==*firstRigidBC) *firstRigidBC=NULL;
	}
	
	// delete any rigid BCs that were not reused
	while(nextBC!=NULL)
	{	prevBC=nextBC;
		nextBC=(BoundaryCondition *)prevBC->GetNextObject();
		delete prevBC;
	}
}

// Debugging aid to count BCs on each step
void ProjectRigidBCsTask::CountBCs(BoundaryCondition **firstBC,BoundaryCondition **lastBC,BoundaryCondition **firstRigidBC)
{
	if(*firstBC==NULL) return;
	
	BoundaryCondition *nextBC=*firstBC;
	int count=0,fullcount=0;
	while(nextBC!=NULL)
	{	if(nextBC==*firstRigidBC)
		{	cout << "# " << count << " on grid, ";
			count=0;
		}
		count++;
		fullcount++;
		if(nextBC==*lastBC)
		{	cout << count << " rigid, ";
			count=0;
		}
		nextBC=(BoundaryCondition *)nextBC->GetNextObject();
	}
	cout << count << " undeleted, total = " << fullcount << endl;
}
