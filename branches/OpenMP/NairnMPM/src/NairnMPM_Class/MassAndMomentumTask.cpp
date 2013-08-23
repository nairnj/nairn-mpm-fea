/********************************************************************************
	MassAndMomentumTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	This first MPM step projects mass and momentum to grid and optionally a few
		other things. It also sets up rigid particle boundary conditions
 
	For each non-rigid particle:
		1. Update dilated volume on the particle (if needed)
		2. Determine crack field
		3. Extrapolate mass and momentum to the grid
		4. If needed extrapolate displacement and unscaled volume
		5. For multimaterial mode, extrapolate mass gradient
 
	For rigid multimaterial particles
		1. Find their current velocity
		2. Other the same (except mass treated different and transport is skipped)
 
	For rigid BC particles
		1. Create the needed boundary condition on the grid
 
	When extrapolation is done, calculate total nodal masses and finish up
		transport tasks.
 
	Finally, adjust momenta for material contact, crack contact, and grid
		bounady conditions
 
	Input Variables
		theElements[], firstCrack..., transportTask...
		mpm[]->vel, pTemperature, pConcentration, volume
		theMaterial[]->GetHeatCapacity()
		
	Output Variables
		mpm[]->ncpos
		mpm[]->vel of rigid particles (only in multimaterial mode)
		mpm[]->volume (dilated volume if needed)
		mpm[]->vfld[]
		cvf[] and their mvf[] are created during initilization
			mvf[]->pk, vk, numberPoints, mass
			mvf[]->disp, volume (if contact might happen)
			mvf[]->volumeGrad (if multimaterial mode)
			cvf[]->norm, numberMaterials
		nd[]->mass
		nd[]->gConcentration, gVolume (for diffusion)
		nd[]->gTemperature, gMpCp (for conduction)
********************************************************************************/

#include "NairnMPM_Class/MassAndMomentumTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/RigidMaterial.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Custom_Tasks/TransportTask.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Exceptions/CommonException.hpp"
#include "Boundary_Conditions/NodalTempBC.hpp"
#include "Boundary_Conditions/NodalConcBC.hpp"
#include "Cracks/CrackNode.hpp"
#include "Nodes/MaterialInterfaceNode.hpp"
#include "Patches/GridPatch.hpp"

#ifdef LOG_PROGRESS
#include "System/ArchiveData.hpp"
#endif

// to ignore crack interactions (only valid if 1 crack or non-interacting cracks)
//#define IGNORE_CRACK_INTERACTIONS

#pragma mark CONSTRUCTORS

MassAndMomentumTask::MassAndMomentumTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Get mass matrix, find dimensionless particle locations,
//	and find grid momenta
void MassAndMomentumTask::Execute(void)
{   
	CommonException *massErr = NULL;
    double fn[maxShapeNodes],xDeriv[maxShapeNodes],yDeriv[maxShapeNodes],zDeriv[maxShapeNodes];
    int nds[maxShapeNodes];
    
    // Set rigid BC contact material velocities first (so loop can be parallel)
	// GetVectorSetting() uses globals and therefore can't be parallel
    if(nmpmsRC>nmpmsNR)
    {   Vector newvel;
        bool hasDir[3];
        for(int p=nmpmsNR;p<nmpmsRC;p++)
        {   MPMBase *mpmptr = mpm[p];
            const RigidMaterial *matID = (RigidMaterial *)theMaterials[mpm[p]->MatID()];
            if(matID->GetVectorSetting(&newvel,hasDir,mtime,&mpmptr->pos))
            {   // change velocity if functions being used, otherwise keep velocity constant
                if(hasDir[0]) mpmptr->vel.x = newvel.x;
                if(hasDir[1]) mpmptr->vel.y = newvel.y;
                if(hasDir[2]) mpmptr->vel.z = newvel.z;
            }
        }
    }
	
	// loop over non-rigid and rigid contact particles - this parallel part changes only particle p
	// mass, momenta, etc are stored on ghost nodes, which are sent to real nodes in next non-parallel loop
    //for(int pn=0;pn<4;pn++)
#pragma omp parallel private(fn,xDeriv,yDeriv,zDeriv,nds)
	{
        // thread for patch pn
		int pn = GetPatchNumber();
        
		// in case 2D planar
        for(int i=0;i<maxShapeNodes;i++) zDeriv[i] = 0.;
        
		try
		{	for(int block=FIRST_NONRIGID;block<=FIRST_RIGID_CONTACT;block++)
			{	MPMBase *mpmptr = patches[pn]->GetFirstBlockPointer(block);
				while(mpmptr!=NULL)
                {   const MaterialBase *matID = theMaterials[mpmptr->MatID()];		// material object for this particle
					int matfld = matID->GetField();									// material velocity field
					
					// get nodes and shape function for material point p
					int i,numnds;
					const ElementBase *elref = theElements[mpmptr->ElemID()];		// element containing this particle
					if(fmobj->multiMaterialMode)
						elref->GetShapeGradients(&numnds,fn,nds,xDeriv,yDeriv,zDeriv,mpmptr);
					else
						elref->GetShapeFunctions(&numnds,fn,nds,mpmptr);
					
					// Add particle property to each node in the element
					short vfld;
					NodalPoint *ndptr;
					for(i=1;i<=numnds;i++)
					{   // get node pointer
                        ndptr = GetNodePointer(pn,nds[i]);
						
						// momentum vector (and allocate velocity field if needed)
						vfld = mpmptr->vfld[i];
						ndptr->AddMassMomentum(mpmptr,vfld,matfld,fn[i],xDeriv[i],yDeriv[i],zDeriv[i],
											   1,block==FIRST_NONRIGID);
					}
					
					// next material point
					mpmptr = (MPMBase *)mpmptr->GetNextObject();
				}
			}
		}
		catch(CommonException err)
        {   if(massErr==NULL)
			{
#pragma omp critical
				massErr = new CommonException(err);
			}
		}
	}
	
	// throw now - only possible error is too many CPDI nodes in 3D
	if(massErr!=NULL) throw *massErr;
    
	// reduction of ghost node forces to real nodes
	int totalPatches = fmobj->GetTotalNumberOfPatches();
	if(totalPatches>1)
	{	for(int pn=0;pn<totalPatches;pn++)
			patches[pn]->MassAndMomentumReduction();
	}
    	
	// undo dynamic velocity, temp, and conc BCs from rigid materials
	UnsetRigidBCs((BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
				  (BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
	UnsetRigidBCs((BoundaryCondition **)&firstTempBC,(BoundaryCondition **)&lastTempBC,
				  (BoundaryCondition **)&firstRigidTempBC,(BoundaryCondition **)&reuseRigidTempBC);
	UnsetRigidBCs((BoundaryCondition **)&firstConcBC,(BoundaryCondition **)&lastConcBC,
				  (BoundaryCondition **)&firstRigidConcBC,(BoundaryCondition **)&reuseRigidConcBC);
	
	// For Rigid BC materials create velocity BC on each node in the element
	for(int p=nmpmsRC;p<nmpms;p++)
	{	MPMBase *mpmptr = mpm[p];										// pointer
		
		int matid0 = mpmptr->MatID();
		const MaterialBase *matID = theMaterials[matid0];				// material object for this particle
		RigidMaterial *rigid=(RigidMaterial *)matID;

		const ElementBase *elref = theElements[mpmptr->ElemID()];		// element containing this particle
		int numnds=elref->NumberNodes();
		
		double rvalue;
		for(int i=1;i<=numnds;i++)
		{   int mi=elref->nodes[i-1];		// 1 based node
			
			// look for setting function in one to three directions
			// GetVectorSetting() returns true if function has set the velocity, otherwise it return FALSE
			bool hasDir[3];
			Vector rvel;
			if(rigid->GetVectorSetting(&rvel,hasDir,mtime,&mpmptr->pos))
			{   // velocity set by 1 to 3 functions as determined by hasDir[i]
				if(hasDir[0])
				{	mpmptr->vel.x = rvel.x;
					SetRigidBCs(mi,matid0,X_DIRECTION,rvel.x,0.,
							(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
							(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
				}
				if(hasDir[1])
				{	mpmptr->vel.y = rvel.y;
					SetRigidBCs(mi,matid0,Y_DIRECTION,rvel.y,0.,
								(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
								(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
				}
				if(hasDir[2])
				{	mpmptr->vel.z = rvel.z;
					SetRigidBCs(mi,matid0,Z_DIRECTION,rvel.z,0.,
								(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
								(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
				}
			}
			else
			{   // velocity set by particle velocity in selected directions
				if(rigid->RigidDirection(X_DIRECTION))
				{	SetRigidBCs(mi,matid0,X_DIRECTION,mpmptr->vel.x,0.,
									(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
									(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
				}
				if(rigid->RigidDirection(Y_DIRECTION))
				{	SetRigidBCs(mi,matid0,Y_DIRECTION,mpmptr->vel.y,0.,
								(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
								(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
				}
				if(rigid->RigidDirection(Z_DIRECTION))
				{	SetRigidBCs(mi,matid0,Z_DIRECTION,mpmptr->vel.z,0.,
								(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
								(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
				}
			}
			
			// temperature
			if(rigid->RigidTemperature())
			{	if(rigid->GetValueSetting(&rvalue,mtime,&mpmptr->pos)) mpmptr->pTemperature=rvalue;
				SetRigidBCs(mi,matid0,TEMP_DIRECTION,mpmptr->pTemperature,0.,
							(BoundaryCondition **)&firstTempBC,(BoundaryCondition **)&lastTempBC,
							(BoundaryCondition **)&firstRigidTempBC,(BoundaryCondition **)&reuseRigidTempBC);
			}
			
			// concentration
			if(rigid->RigidConcentration())
			{	if(rigid->GetValueSetting(&rvalue,mtime,&mpmptr->pos)) mpmptr->pConcentration=rvalue;
				SetRigidBCs(mi,matid0,CONC_DIRECTION,mpmptr->pConcentration,0.,
							(BoundaryCondition **)&firstConcBC,(BoundaryCondition **)&lastConcBC,
							(BoundaryCondition **)&firstRigidConcBC,(BoundaryCondition **)&reuseRigidConcBC);
			}
		}
	}
	
	// if any left over rigid BCs, delete them now
	RemoveRigidBCs((BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,(BoundaryCondition **)&firstRigidVelocityBC);
	RemoveRigidBCs((BoundaryCondition **)&firstTempBC,(BoundaryCondition **)&lastTempBC,(BoundaryCondition **)&firstRigidTempBC);
	RemoveRigidBCs((BoundaryCondition **)&firstConcBC,(BoundaryCondition **)&lastConcBC,(BoundaryCondition **)&firstRigidConcBC);
	
#ifdef COMBINE_RIGID_MATERIALS
	bool combineRigid = firstCrack!=NULL && fmobj->multiMaterialMode && fmobj->hasRigidContactParticles;
#endif
	
	// Post mass and momentum extrapolation calculations on nodes
#pragma omp parallel
	{
		// variables for each thread
		CrackNode *firstCrackNode=NULL,*lastCrackNode=NULL;
		MaterialInterfaceNode *firstInterfaceNode=NULL,*lastInterfaceNode=NULL;
		
		// Each pass in this loop should be independent
#pragma omp for nowait
		for(int i=1;i<=nnodes;i++)
		{	// node reference
			NodalPoint *ndptr = nd[i];
			
			try
            {
#ifdef COMBINE_RIGID_MATERIALS
                // combine rigid fields if necessary
                if(combineRigid)
                    ndptr->CopyRigidParticleField();
#endif
				// Get total nodal masses and count materials if multimaterial mode
				ndptr->CalcTotalMassAndCount();

				// multimaterial contact
				if(fmobj->multiMaterialMode)
					ndptr->MaterialContactOnNode(timestep,MASS_MOMENTUM_CALL,&firstInterfaceNode,&lastInterfaceNode);
				
				// crack contact
				if(firstCrack!=NULL)
					ndptr->CrackContact(FALSE,0.,&firstCrackNode,&lastCrackNode);
				
				// get transport values on nodes
				TransportTask *nextTransport=transportTasks;
				while(nextTransport!=NULL)
					nextTransport = nextTransport->GetNodalValue(ndptr);
			}
			catch(CommonException err)
			{	if(massErr==NULL)
				{
#pragma omp critical
					massErr = new CommonException(err);
				}
			}
		}

#pragma omp critical
		{
			// link up crack nodes
			if(lastCrackNode != NULL)
			{	if(CrackNode::currentCNode != NULL)
					firstCrackNode->SetPrevBC(CrackNode::currentCNode);
				CrackNode::currentCNode = lastCrackNode;
			}
			
			// link up interface nodes
			if(lastInterfaceNode != NULL)
			{	if(MaterialInterfaceNode::currentIntNode != NULL)
					firstInterfaceNode->SetPrevBC(MaterialInterfaceNode::currentIntNode);
				MaterialInterfaceNode::currentIntNode = lastInterfaceNode;
			}
		}
	}
	
	// throw any errors
	if(massErr!=NULL) throw *massErr;
    
	// Impose transport BCs and extrapolate gradients to the particles
	TransportTask *nextTransport=transportTasks;
	while(nextTransport!=NULL)
    {   nextTransport->ImposeValueBCs(mtime);
		nextTransport = nextTransport->GetGradients(mtime);
	}
	
	// used to call class methods for material contact and crack contact here
	// Impose velocity BCs
	NodalVelBC::GridMomentumConditions(TRUE);

}

// Set boundary conditions determined by moving rigid paticles
void MassAndMomentumTask::SetRigidBCs(int mi,int matid0,int type,double value,double angle,BoundaryCondition **firstBC,
						   BoundaryCondition **lastBC,BoundaryCondition **firstRigidBC,BoundaryCondition **reuseRigidBC)
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
            {   int newType = type == Z_DIRECTION ? Z_DIRECTION_INPUT : type ;
				newBC=(BoundaryCondition *)(new NodalVelBC(mi,newType,CONSTANT_VALUE,value,(double)0.,(double)0.,(double)0.));
				if(newBC==NULL) throw CommonException("Memory error allocating rigid particle boundary condition.",
													  "NairnMPM::SetRigidBCs");
			}
			break;
			
		case TEMP_DIRECTION:
			if(*reuseRigidBC!=NULL)
				newBC=(*reuseRigidBC)->SetRigidProperties(mi,type,CONSTANT_VALUE,value);
			else
			{	newBC=(BoundaryCondition *)(new NodalTempBC(mi,CONSTANT_VALUE,value,(double)0.));
				if(newBC==NULL) throw CommonException("Memory error allocating rigid particle boundary condition.",
													  "NairnMPM::SetRigidBCs");
			}
			break;
			
		case CONC_DIRECTION:
			if(*reuseRigidBC!=NULL)
				newBC=(*reuseRigidBC)->SetRigidProperties(mi,type,CONSTANT_VALUE,value);
			else
			{	newBC=(BoundaryCondition *)(new NodalConcBC(mi,CONSTANT_VALUE,value,(double)0.));
				if(newBC==NULL) throw CommonException("Memory error allocating rigid particle boundary condition.",
													  "NairnMPM::SetRigidBCs");
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
// that are no longer needed
void MassAndMomentumTask::RemoveRigidBCs(BoundaryCondition **firstBC,BoundaryCondition **lastBC,BoundaryCondition **firstRigidBC)
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

// Unset nodal dof for rigid BCs so can try to reuse
// them without needing new memory allocations
void MassAndMomentumTask::UnsetRigidBCs(BoundaryCondition **firstBC,BoundaryCondition **lastBC,
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

// Debugging aid to count BCs on each step
void MassAndMomentumTask::CountBCs(BoundaryCondition **firstBC,BoundaryCondition **lastBC,BoundaryCondition **firstRigidBC)
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


