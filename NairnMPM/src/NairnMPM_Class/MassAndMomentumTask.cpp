/********************************************************************************
	MassAndMomentumTask.cpp
	NairnMPM

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
		theMaterial[]->GetHeatCapacity() (for conduction)
		
	Output Variables
		mpm[]->ncpos
		mpm[]->vel of rigid particles (only in multimaterial mode)
		mpm[]->volume (dilated volume if needed)
		mpm[]->vfld[]
		Allocate cvf[] and their mvf[] as needed
			mvf[]->pk, vk, numberPoints, mass
			mvf[]->disp, unscaledVolume (if contact might happen)
			mvf[]->massGrad (if multimaterial mode)
			cvf[]->norm, numberMaterials
		nd[]->mass
		nd[]->gConcentration, gVolume (for diffusion)
		nd[]->gTemperature, gRhoVCp (for conduction)
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

// to ignore crack interactions (only valid if 1 crack or non-interacting cracks)
//#define IGNORE_CRACK_INTERACTIONS

// uncomment to project rigid velocity fields to all crack velocity fields
// Only does anything when have cracks, in multimaterial mode, and has rigid contact particles
//#define COMBINE_RIGID_MATERIALS

#pragma mark CONSTRUCTORS

MassAndMomentumTask::MassAndMomentumTask(const char *name) : MPMTask(name)
{
	// zero thisfunction in case in 2D analysis
	for(int i=0;i<MaxShapeNds;i++) zDeriv[i]=0.;
	
	// crack locations defaults (in case no cracks)
	cfld[0].loc=NO_CRACK;
	cfld[1].loc=NO_CRACK;
}

#pragma mark REQUIRED METHODS

// Get mass matrix, find dimensionless particle locations,
//	and find grid momenta
void MassAndMomentumTask::Execute(void)
{
#ifdef _PROFILE_TASKS_
	double beginTime=fmobj->CPUTime();
#endif
	
	int i,p,mi,iel,matfld,numnds,nds[MaxShapeNds];
	MaterialBase *matID;
	double mp,fn[MaxShapeNds],xDeriv[MaxShapeNds],yDeriv[MaxShapeNds];
    CrackHeader *nextCrack;
    short vfld;
	TransportTask *nextTransport;
	MPMBase *mpmptr;
	NodalPoint *ndptr;
	
	// undo dynamic velocity, temp, and conc BCs from rigid materials
	UnsetRigidBCs((BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
				  (BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
	UnsetRigidBCs((BoundaryCondition **)&firstTempBC,(BoundaryCondition **)&lastTempBC,
				  (BoundaryCondition **)&firstRigidTempBC,(BoundaryCondition **)&reuseRigidTempBC);
	UnsetRigidBCs((BoundaryCondition **)&firstConcBC,(BoundaryCondition **)&lastConcBC,
				  (BoundaryCondition **)&firstRigidConcBC,(BoundaryCondition **)&reuseRigidConcBC);
	
	// loop over particles
    for(p=0;p<nmpms;p++)
	{	mpmptr=mpm[p];							// pointer
		iel=mpmptr->ElemID();					// element containing this particle
		matID=theMaterials[mpmptr->MatID()];	// material object for this particle
		
		// normal materials
		if(!matID->RigidBC())
		{	mp=mpmptr->mp;					// material point mass in g
			matfld=matID->GetField();		// material velocity field
			
			// get nodes and shape function for material point p
			if(fmobj->multiMaterialMode)
			{	theElements[iel]->GetShapeFunctionsAndGradients(&numnds,fn,nds,&mpmptr->pos,mpmptr->GetNcpos(),xDeriv,yDeriv,zDeriv);
				
				// for particles that are multimaterial rigid materials, set their velocity
				if(matID->Rigid())
				{	Vector newvel;
					bool hasDir[3];
					if(((RigidMaterial *)matID)->GetVectorSetting(&newvel,hasDir,mtime,&mpmptr->pos))
					{	if(hasDir[0]) mpmptr->vel.x = newvel.x;
						if(hasDir[1]) mpmptr->vel.y = newvel.y;
						if(hasDir[2]) mpmptr->vel.z = newvel.z;
					}
				}
			}
			else
				theElements[iel]->GetShapeFunctions(&numnds,fn,nds,&mpmptr->pos,mpmptr->GetNcpos());
			
			// get deformed particle volume if it will be needed (for transport tasks)
			if(fmobj->volumeExtrap) mpmptr->SetDilatedVolume();
			
			// Add particle property to each node in the element
			for(i=1;i<=numnds;i++)
			{	ndptr=nd[nds[i]];				// get pointer
				
				// Look for crack crossing and save until later
				if(firstCrack!=NULL)
				{	int cfound=0;
					Vector norm;
					cfld[0].loc=NO_CRACK;			// NO_CRACK, ABOVE_CRACK, or BELOW_CRACK
					cfld[1].loc=NO_CRACK;
					nextCrack=firstCrack;
					while(nextCrack!=NULL)
					{   vfld=nextCrack->CrackCross(mpmptr->pos.x,mpmptr->pos.y,ndptr->x,ndptr->y,&norm);
						if(vfld!=NO_CRACK)
						{	cfld[cfound].loc=vfld;
							cfld[cfound].norm=norm;
#ifdef IGNORE_CRACK_INTERACTIONS
							cfld[cfound].crackNum=1;	// appears to always be same crack, and stop when found one
							break;
#else
							cfld[cfound].crackNum=nextCrack->GetNumber();
							cfound++;
							if(cfound>1) break;			// stop if found two, if there are more then two, physics will be off
#endif
						}
						nextCrack=(CrackHeader *)nextCrack->GetNextObject();
					}
				}
				
				// momentum vector (and allocate velocity field if needed)
				vfld=ndptr->AddMomentumTask1(matfld,cfld,fn[i]*mp,&mpmptr->vel);
				mpmptr->vfld[i]=vfld;
				
				// crack contact calculations
				contact.AddDisplacementVolumeTask1(vfld,matfld,ndptr,mpmptr,fn[i]);
				
				// material contact calculations
				if(fmobj->multiMaterialMode)
					ndptr->AddMassGradient(vfld,matfld,mp,xDeriv[i],yDeriv[i],zDeriv[i]);
				
				// more for non-rigid contact materials
				if(!matID->Rigid())
				{	// add to lumped mass matrix
					ndptr->AddMass(vfld,matfld,mp*fn[i]);
					
					// transport calculations
					nextTransport=transportTasks;
					while(nextTransport!=NULL)
						nextTransport=nextTransport->Task1Extrapolation(ndptr,mpmptr,fn[i]);
				}
				else
				{	// for rigid particles, let the crack velocity field know
					ndptr->AddMassTask1(vfld,matfld);
				}
			}
		}
		
		// For Rigid BC materials create velocity BC on each node in the element
		else
		{	numnds=theElements[iel]->NumberNodes();
			double rvalue;
			RigidMaterial *rigid=(RigidMaterial *)matID;
			for(i=1;i<=numnds;i++)
			{   mi=theElements[iel]->nodes[i-1];		// 1 based node
				
				// look for setting function in one to three directions
				bool hasDir[3];
				Vector rvel;
				if(rigid->GetVectorSetting(&rvel,hasDir,mtime,&mpmptr->pos))
				{	if(hasDir[0])
					{	mpmptr->vel.x = rvel.x;
						SetRigidBCs(mi,X_DIRECTION,rvel.x,0.,
								(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
								(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
					}
					if(hasDir[1])
					{	mpmptr->vel.y = rvel.y;
						SetRigidBCs(mi,Y_DIRECTION,rvel.y,0.,
									(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
									(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
					}
					if(hasDir[2])
					{	mpmptr->vel.z = rvel.z;
						SetRigidBCs(mi,Z_DIRECTION,rvel.z,0.,
									(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
									(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
					}
				}
				else
				{	// check skewed, x or y direction velocities
					if(rigid->RigidDirection(X_DIRECTION))
					{	SetRigidBCs(mi,X_DIRECTION,mpmptr->vel.x,0.,
										(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
										(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
					}
					if(rigid->RigidDirection(Y_DIRECTION))
					{	SetRigidBCs(mi,Y_DIRECTION,mpmptr->vel.y,0.,
									(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
									(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
					}
					if(rigid->RigidDirection(Z_DIRECTION))
					{	SetRigidBCs(mi,Z_DIRECTION,mpmptr->vel.z,0.,
									(BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,
									(BoundaryCondition **)&firstRigidVelocityBC,(BoundaryCondition **)&reuseRigidVelocityBC);
					}
				}
				
				// temperature
				if(rigid->RigidTemperature())
				{	if(rigid->GetSetting(&rvalue,mtime,&mpmptr->pos)) mpmptr->pTemperature=rvalue;
					SetRigidBCs(mi,TEMP_DIRECTION,mpmptr->pTemperature,0.,
								(BoundaryCondition **)&firstTempBC,(BoundaryCondition **)&lastTempBC,
								(BoundaryCondition **)&firstRigidTempBC,(BoundaryCondition **)&reuseRigidTempBC);
				}
				
				// concentration
				if(rigid->RigidConcentration())
				{	if(rigid->GetSetting(&rvalue,mtime,&mpmptr->pos)) mpmptr->pConcentration=rvalue;
					SetRigidBCs(mi,CONC_DIRECTION,mpmptr->pConcentration,0.,
								(BoundaryCondition **)&firstConcBC,(BoundaryCondition **)&lastConcBC,
								(BoundaryCondition **)&firstRigidConcBC,(BoundaryCondition **)&reuseRigidConcBC);
				}
			}
		}
    }
	
	// if any left over rigid BCs, delete them now
	RemoveRigidBCs((BoundaryCondition **)&firstVelocityBC,(BoundaryCondition **)&lastVelocityBC,(BoundaryCondition **)&firstRigidVelocityBC);
	RemoveRigidBCs((BoundaryCondition **)&firstTempBC,(BoundaryCondition **)&lastTempBC,(BoundaryCondition **)&firstRigidTempBC);
	RemoveRigidBCs((BoundaryCondition **)&firstConcBC,(BoundaryCondition **)&lastConcBC,(BoundaryCondition **)&firstRigidConcBC);
	
	// Get total nodal masses and count materials if multimaterial mode
	for(i=1;i<=nnodes;i++)
		nd[i]->CalcTotalMassAndCount();
	
#ifdef COMBINE_RIGID_MATERIALS
	// combine rigid fields if necessary
	if(firstCrack!=NULL && fmobj->multiMaterialMode && fmobj->hasRigidContactParticles)
		NodalPoint::CombineRigidMaterials();
#endif
	
	// Find values and gradients for transport tasks
	nextTransport=transportTasks;
	while(nextTransport!=NULL)
		nextTransport=nextTransport->GetValuesAndGradients(mtime);
	
	// Adjust momenta for multimaterial contact
	NodalPoint::MaterialContact(fmobj->multiMaterialMode,FALSE,timestep);
	
	// Adjust momenta for crack contact
	CrackHeader::ContactConditions(TRUE);
	
	// Impose velocity BCs
	// NOTE: Switched order of contact and BCs (8/12/2009)
	NodalVelBC::GridMomentumConditions(TRUE);
	
#ifdef _PROFILE_TASKS_
	totalTaskTime+=fmobj->CPUTime()-beginTime;
#endif
}	

// Set boundary conditions determined by moving rigid paticles
void MassAndMomentumTask::SetRigidBCs(int mi,int type,double value,double angle,BoundaryCondition **firstBC,
						   BoundaryCondition **lastBC,BoundaryCondition **firstRigidBC,BoundaryCondition **reuseRigidBC)
{
	BoundaryCondition *newBC=NULL;
	NodalVelBC *velBC;
	
	// check if already set in that direction by actual BC or by previous rigid BC
	// New rigid BC's can only be on free directions
	if(nd[mi]->fixedDirection&type) return;
	
	// create new boundary conditions
	switch(type)
	{	case SKEW_DIRECTION:
			if(nd[mi]->fixedDirection&(X_DIRECTION+Y_DIRECTION)) return;
			if(*reuseRigidBC!=NULL)
				velBC=(NodalVelBC *)((*reuseRigidBC)->SetRigidProperties(mi,type,CONSTANT_VALUE,value));
			else
			{	velBC=new NodalVelBC(mi,type,CONSTANT_VALUE,value,(double)0.);
				if(velBC==NULL) throw CommonException("Memory error allocating rigid particle boundary condition.",
													  "NairnMPM::SetRigidBCs");
			}
			velBC->SetSkewAngle(angle);
			newBC=(BoundaryCondition *)velBC;
			break;
			
		case X_DIRECTION:
		case Y_DIRECTION:
		case Z_DIRECTION:
			if(*reuseRigidBC!=NULL)
				newBC=(*reuseRigidBC)->SetRigidProperties(mi,type,CONSTANT_VALUE,value);
			else
			{	newBC=(BoundaryCondition *)(new NodalVelBC(mi,type,CONSTANT_VALUE,value,(double)0.));
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


