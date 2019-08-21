/********************************************************************************
	CrackVelocityFieldSingle.cpp
	NairnMPM

	Created by John Nairn on 21 August 2009.
	Copyright (c) 2009 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Nodes/CrackVelocityFieldSingle.hpp"
#include "Exceptions/CommonException.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Global_Quantities/BodyForce.hpp"

#pragma mark INITIALIZATION

// constructor
CrackVelocityFieldSingle::CrackVelocityFieldSingle(int num,short theLoc,int cnum) : CrackVelocityField(num,theLoc,cnum)
{
}

// Destructor
CrackVelocityFieldSingle::~CrackVelocityFieldSingle()
{	// delete velocity field if present
	if(mvf[0]!=NULL) delete mvf[0];
}

// zero the one material if it was active
void CrackVelocityFieldSingle::ZeroMatFields(void)
{	mvf[0]->Zero();
}

#pragma mark TASK 1 METHODS

// Make copy of momenta, to be restored after strain update and forces and return mass
// Note: never have materials that ignore cracks when this class is active
double CrackVelocityFieldSingle::GetTotalMassAndCount(bool &hasMaterialContact)
{	if(mvf[0]->numberPoints<1) return 0.;
	return mvf[0]->GetTotalMassAndCount();
}

#pragma mark TASK 3 METHODS

// Add to force spread out over the materials so each has same extra accelerations = f/M_i
// Only called by AddTractionForce()
// Only adds force to fields that see cracks (not relevant in single material mode)
void CrackVelocityFieldSingle::AddFtotSpreadTask3(Vector *f)
{	if(mvf[0]->numberPoints>0)
		mvf[0]->AddFtot(f);
}

// Add gravity and body force at a node in g mm/sec^2
void CrackVelocityFieldSingle::AddGravityAndBodyForceTask3(Vector *gridBodyForce,double gridAlpha,double gridForceAlpha)
{	if(mvf[0]->numberPoints>0)
		mvf[0]->AddGravityAndBodyForceTask3(gridBodyForce,gridAlpha,gridForceAlpha);
}

// Restore momenta just prior to momentum update and to setting forces
// for velocity BCs
void CrackVelocityFieldSingle::RestoreMomenta(void)
{	if(mvf[0]->numberPoints>0)
	{	mvf[0]->RestoreMomenta();
	}
}

#pragma mark TASK 4 METHODS

// update momenta for this MPM step
//  pk(i+1) = pk(i) + ftot * dt
void CrackVelocityFieldSingle::UpdateMomentum(double timestep)
{	if(mvf[0]->numberPoints>0)
        mvf[0]->UpdateMomentum(timestep);
}

#pragma mark TASK 5 XPIC METHODS

// Support XPIC calculations
void CrackVelocityFieldSingle::XPICSupport(int xpicCalculation,int xpicOption,NodalPoint *real,double timestep,int m,int k,double vsign)
{	if(xpicCalculation==COPY_VSTARNEXT)
	{	xpicOption = fieldNum;			// vfld
		m = 0;							// matfld
	}
	mvf[0]->XPICSupport(xpicCalculation,xpicOption,real,timestep,m,k,vsign);
}

#pragma mark TASK 6 METHODS

// zero momentum and displacement at a node for new calculations
void CrackVelocityFieldSingle::RezeroNodeTask6(double delTime)
{	if(mvf[0]->numberPoints>0)
		mvf[0]->RezeroNodeTask6();
}

#pragma mark VELOCITY METHODS

// Various calculations of grid values
void CrackVelocityFieldSingle::GridValueCalculation(int calcOption)
{	if(mvf[0]->numberPoints>0)
		mvf[0]->GridValueCalculation(calcOption);
}

#pragma mark BOUNDARY CONDITIONS

// zero one component of moment and velocity on one material field
void CrackVelocityFieldSingle::ZeroVelocityBC(Vector *norm,int passType,double deltime,Vector *freaction)
{	if(mvf[0]->numberPoints>0)
        mvf[0]->ZeroVelocityBC(norm,passType,deltime,freaction);
}

// add one component of momentum and velocity from BCs to one material field
void CrackVelocityFieldSingle::AddVelocityBC(Vector *norm,double vel,int passType,double deltime,Vector *freaction)
{	if(mvf[0]->numberPoints>0)
        mvf[0]->AddVelocityBC(norm,vel,passType,deltime,freaction);
}

// Reflect one component of velocity and momentum from a node
void CrackVelocityFieldSingle::ReflectVelocityBC(Vector *norm,CrackVelocityField *rcvf,double vel0,double reflectRatio,
												 int passType,double deltime,Vector *freaction)
{	if(mvf[0]->numberPoints>0)
	{	MatVelocityField **rmvf = rcvf->GetMaterialVelocityFields();
		if(rmvf[0]->numberPoints>0)
		{	double rvel;
			if(passType==UPDATE_GRID_STRAINS_CALL)
				rvel = vel0 + reflectRatio*(vel0 - DotVectors(norm,rmvf[0]->vk));
			else
				rvel = vel0 + reflectRatio*(vel0 - DotVectors(norm,&rmvf[0]->pk)/rmvf[0]->mass);
			mvf[0]->AddVelocityBC(norm,rvel,passType,deltime,freaction);
		}
	}
}

#pragma mark CONTACT ACCESSORS

// get volume when only a single material (overridden when might be more)
double CrackVelocityFieldSingle::GetContactVolumeNonrigid(bool requireCracks) const
{	return mvf[0]->numberPoints>0 ? mvf[0]->GetContactVolume() : 0. ;
}

// get center of mass displacement (actually sum of displacement*mass so displacement is vector/total mass)
// If on symmetry plane, that component will be zeroed out
// Note: in single material mode, following two are identical
Vector CrackVelocityFieldSingle::GetCMDisplacement(NodalPoint *np,bool requireCracks,bool useDisps) const
{   Vector dk = *(mvf[0]->GetContactDispPtr(useDisps));
	AdjustForSymmetry(np,&dk,false);
	return dk;
}

// get force for the one material fields in this crack velocity field
Vector CrackVelocityFieldSingle::GetCMatFtot(void)
{   if(mvf[0]->numberPoints>0)
		return mvf[0]->GetFtot();
	else
	{	Vector fk;
		ZeroVector(&fk);
		return fk;
	}
}

#pragma mark ACCESSORS

// total mass all velocity fields, but only one in single material mode
// Note: this class is not active when any materials might ignore cracks
// Note: in single material mode, following two are identical
double CrackVelocityFieldSingle::GetTotalMass(bool requireCracks) const
{	return mvf[0]->numberPoints>0 ? mvf[0]->mass : 0. ;
}

// total mass and kinetric energy all velocity fields (rigid particles not counted)
// in g-mm2/sec^2 = nanoJ
void CrackVelocityFieldSingle::AddKineticEnergyAndMass(double &kineticEnergy,double &totalMass)
{	if(mvf[0]->numberPoints>0 && !mvf[0]->IsRigidField())
	{	totalMass += mvf[0]->mass;
		double magp = DotVectors(&mvf[0]->pk,&mvf[0]->pk);
		kineticEnergy += 0.5*magp/mvf[0]->mass;
	}
}

// get center of mass momentum for all material fields in this crack velocity field
// if totalFtot!=NULL, set it to nodal fi = mi ai
Vector CrackVelocityFieldSingle::GetCMatMomentum(bool &hasParticles,double *fieldMass,Vector *totalFtot,bool useVelocity) const
{	hasParticles = mvf[0]->numberPoints>0 ;
	*fieldMass = mvf[0]->mass;
	if(totalFtot!=NULL) *totalFtot=mvf[0]->GetFtot();
	if(useVelocity)
	{	Vector ptot = SetScaledVector(&mvf[0]->vk[0],mvf[0]->mass);
		return ptot;
	}
	else
		return mvf[0]->pk ;
}

/* in response to crack contact, change moment by changing velocity of all materials
	the same amount
 
   Change velocity by dP/M, where M is total mass
   Material i velocity becomes vi = pi/mi + dP/M
   Material i momentum change is mi vi = pi + mi dP/M
*/
void CrackVelocityFieldSingle::ChangeCrackMomentum(Vector *delP,int callType,double deltime)
{	if(mvf[0]->numberPoints>0)
		mvf[0]->ChangeMatMomentum(delP,callType,deltime);
}

// copy all material velocity fields for boundary conditions methods, returning new offset into the save array
int CrackVelocityFieldSingle::CopyFieldMomenta(Vector *holdPk,int offset)
{	if(mvf[0]->numberPoints>0)
	{	holdPk[offset]=mvf[0]->pk;
		offset++;
	}
	return offset;
}

#if ADJUST_COPIED_PK == 1
// set symmetry plane momenta to zero
void CrackVelocityFieldSingle::AdjustForSymmetryBC(NodalPoint *ndptr)
{	if(mvf[0]->numberPoints>0)
    {	mvf[0]->AdjustForSymmetryBC(ndptr->fixedDirection);
	}
}
#endif

// paste all material velocity fields back for boundary conditions methods, returning new offset into the saved array
int CrackVelocityFieldSingle::PasteFieldMomenta(Vector *holdPk,int offset)
{	if(mvf[0]->numberPoints>0)
	{	mvf[0]->pk=holdPk[offset];
		offset++;
	}
	return offset;
}

// for debugging
void CrackVelocityFieldSingle::Describe(void) const
{	CrackVelocityField::Describe();
	cout << "#     single material" << endl;
	mvf[0]->Describe(0);
}
	
