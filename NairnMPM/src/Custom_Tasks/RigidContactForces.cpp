/******************************************************************************** 
	RigidContactForces.cpp
	NairnMPM

	Created by John Nairn on 7/29/10.
	Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Custom_Tasks/RigidContactForces.hpp"
#include "System/ArchiveData.hpp"
#include "Nodes/NodalPoint.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"

#pragma mark INITIALIZE

// Constructors
RigidContactForces::RigidContactForces() { }

// Return name of this task
const char *RigidContactForces::TaskName(void) { return "Extrapolate contact forces to rigid particles"; }

#pragma mark GENERIC TASK METHODS

// called once at start of MPM analysis - initialize and print info
CustomTask *RigidContactForces::Initialize(void)
{
    cout << "Extrapolating contact forces to rigid particles" << endl;
    return nextTask;
}

// called when MPM step is getting ready to do custom tasks
// return nextTask when done
CustomTask *RigidContactForces::PrepareForStep(bool &doNodalExtraps)
{	if((getForcesThisStep=archiver->WillArchive()))
	{	// only need to extrapolate forces when about to archive
		doNodalExtraps=TRUE;
		ZeroVector(&force);
	}
	return nextTask;
}

#pragma mark TASK EXTRAPOLATION METHODS

// used to extapolate nodal values to particles during extrapolations
CustomTask *RigidContactForces::ParticleCalculation(NodalPoint *ndmi,MPMBase *mpnt,short vfld,int matfld,
											double fn,double xDeriv,double yDeriv,double zDeriv,short isRigid)
{
	if(!getForcesThisStep || !isRigid) return  nextTask;
	
	// extrapolate forces stored in this rigid materials velocity field
	Vector ftot=ndmi->GetContactForce(vfld,matfld);
	AddScaledVector(&force,&ftot,fn);
	
	return nextTask;
}

// add particle data to some calculation
CustomTask *RigidContactForces::ParticleExtrapolation(MPMBase *mpnt,short isRigid)
{
	if(!getForcesThisStep || !isRigid) return nextTask;
		
	// copy to stress as force in units of microN per rho of this material
	// rho of rigid particles is always 1000
	// Must multiply by -1/timestep too to get force, since only momenta were tracked by contact code
	Tensor *sp=mpnt->GetStressTensor();
	double scale=-0.001/timestep;
	sp->xx=scale*force.x;
	sp->yy=scale*force.y;
	sp->zz=scale*force.z;
	return nextTask;
}




