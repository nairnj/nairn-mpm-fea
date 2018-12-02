/********************************************************************************
    MatVelocityField.cpp
    nairn-mpm-fea
    
    Created by John Nairn on 3 April 2009.
    Copyright (c) 2009 John A. Nairn, All rights reserved.
 
	Special case for rigid material
		pk will be sum of Vp*vel
		disp will be sum Vp*disp
		mass grad will be volume gradient
		unscaled volume if the volume
	
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/MatVelocityField.hpp"
#include "Exceptions/CommonException.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Global_Quantities/BodyForce.hpp"

// class statics
int MatVelocityField::pkCopy=0;

#pragma mark INITIALIZATION

// Constructors
// throws std::bad_alloc
MatVelocityField::MatVelocityField(int setFlags)
{	if(fmobj->multiMaterialMode)
		volumeGrad=new Vector;
	else
		volumeGrad=NULL;
    
	// New one vector to copy momenta
	xpic = new  Vector[1];
	pkCopy = 0;
	
	SetRigidField(false);
	Zero();
	flags = setFlags;
}

// Destructor
MatVelocityField::~MatVelocityField()
{	if(volumeGrad!=NULL)
		delete volumeGrad;
	delete [] xpic;
}

// zero data at start of time step
void MatVelocityField::Zero(void)
{	mass=0.;
	ZeroVector(&pk);
	if(!IsFixedRigidField())
    {   // rigid fields are tracking contact force on ftot
        // it is cummulative and therefore not zeroed here
		// but RigidBlock fields do clear it
		ZeroVector(&ftot);
	}
	ZeroVector(&vk);
	ZeroVector(&disp);
	if(volumeGrad!=NULL) ZeroVector(volumeGrad);
	numberPoints=0;
	volume=0.;
	ZeroVector(&xpic[pkCopy]);
}

#pragma mark METHODS

// add to momentum in task 1. Save velocity for the first point in case only
// has one point (for more accurate velocity calculation)
// Momentum is g-mm/sec
void MatVelocityField::AddMomentumTask1(Vector *addPk,Vector *vel,int numPts)
{	AddVector(&pk,addPk);
	numberPoints += numPts;
}

// copy mass and momentum from ghost to real node
void MatVelocityField::CopyMassAndMomentum(NodalPoint *real,int vfld,int matfld)
{
	// skip if none
	if(numberPoints==0) return;
	
	// momentum
	real->AddMomentumTask1((short)vfld,matfld,1.,&pk,numberPoints);
	
	// if cracks and/or multimaterial mode
	if(firstCrack!=NULL || fmobj->multiMaterialMode)
	{	real->AddDisplacement((short)vfld,matfld,1.,&disp);
		real->AddVolume((short)vfld,matfld,volume);
        
        // volume gradient
        if(fmobj->multiMaterialMode)
            real->CopyVolumeGradient((short)vfld,matfld,volumeGrad);
	}
	
	// if rigid particle in multimaterial mode
	if(!IsRigidField())
	{	// mass
		real->AddMass((short)vfld,matfld,mass);
	}
	else
	{	// mass to count particles
		real->AddMassTask1((short)vfld,matfld,mass,numberPoints);
	}
}

// copy mass and momentum from ghost to real node
void MatVelocityField::CopyMassAndMomentumLast(NodalPoint *real,int vfld,int matfld)
{
	// skip if none or rigid
	if(numberPoints==0 || IsRigidField()) return;
	
	// momentum
	real->AddMomentumTask6((short)vfld,matfld,1.,&pk);
	
	// if cracks and multimaterial mode
	if(firstCrack!=NULL || fmobj->multiMaterialMode)
	{	real->AddDisplacement((short)vfld,matfld,1.,&disp);
		real->AddVolume((short)vfld,matfld,volume);
        
        // volume gradient
        if(fmobj->multiMaterialMode)
            real->CopyVolumeGradient((short)vfld,matfld,volumeGrad);
	}
}

// Copy grid forces from this ghost node to the real node (nonrigid only)
void MatVelocityField::CopyGridForces(NodalPoint *real,int vfld,int matfld)
{	
	// skip if none
	if(numberPoints==0 || IsRigidField()) return;
	
	// Ftot
	real->AddFtotTask3(vfld,matfld,&ftot);
}

// Restore momenta
void MatVelocityField::RestoreMomenta(void)
{
	// paste the extrapolated momenta back and zero storage location
	pk = xpic[pkCopy];
}

// in response to contact, change the momentum
// for a single point, calculate the velocity
// if postUpdate is true, then adjust ftot as well to keep in sync with momentum change
// Momentum is g-mm/sec
void MatVelocityField::ChangeMatMomentum(Vector *delP,int callType,double deltime)
{	AddVector(&pk,delP);
	if(callType==UPDATE_MOMENTUM_CALL)
	{	AddScaledVector(&ftot, delP, 1./deltime);
	}
	else if(callType==MASS_MOMENTUM_CALL)
	{	xpic[pkCopy] = pk;
	}
}

// in response to contact, and only for rigid materials, add contact force to ftot
// Contact force is g-mm/sec^2 = micro N (because N is kg-m/sec^2)
// Here just sum momenta applied to non-rigid materials. To get forces, this result
//	is multiplied by -1/timestep
void MatVelocityField::AddContactForce(Vector *delP)
{	ftot.x += delP->x;
	ftot.y += delP->y;
	ftot.z += delP->z;
}

// Calculate velocity at a node from current momentum and mass matrix in all velocity fields
// Velocity is mm/sec
// Note: do not call for rigid fields in multimaterial mode
void MatVelocityField::CalcVelocityForStrainUpdate(void)
{	// exit if nothing there
	if(numberPoints==0 || mass==0.) return;
	
	// get velocity
	CopyScaleVector(&vk,&pk,1./mass);
}

// Add grid damping force at a node in g mm/sec^2
void MatVelocityField::AddGravityAndBodyForceTask3(Vector *gridBodyForce)
{	ftot.x += mass*gridBodyForce->x;
	ftot.y += mass*gridBodyForce->y;
	ftot.z += mass*gridBodyForce->z;
}

// internal force - add or scale and add
void MatVelocityField::AddFtot(Vector *f)
{	ftot.x += f->x;
	ftot.y += f->y;
	ftot.z += f->z;
}
void MatVelocityField::AddFtotScaled(Vector *f,double scaled)
{	ftot.x += f->x*scaled;
	ftot.y += f->y*scaled;
	ftot.z += f->z*scaled;
}

// Update momentum for this MPM step
//  pk(i+1) = pk(i) + ftot * dt
void MatVelocityField::UpdateMomentum(double timestep)
{	// update momenta
    AddScaledVector(&pk,&ftot,timestep);
}

// on particle updates, increment nodal velocity and acceleration and others as needed
// fi is shape function
void MatVelocityField::IncrementNodalVelAcc(double fi,GridToParticleExtrap *gp) const
{
    double mnode = fi/mass;					// Ni/mass
	
	// Summing S v+
	gp->vgpnp1.x += pk.x*mnode;				// velocity += p/mass
	gp->vgpnp1.y += pk.y*mnode;
	gp->vgpnp1.z += pk.z*mnode;
	
	// Summing S a
	gp->acc->x += ftot.x*mnode;				// acceleration += f/mass
	gp->acc->y += ftot.y*mnode;
	gp->acc->z += ftot.z*mnode;
}

#pragma mark ACCESSORS

// for debugging
void MatVelocityField::Describe(int fldnum) const
{
	cout << "#      " << fldnum << ". Material Field: n="<<  numberPoints << " mass=" << mass ;
    cout << " p=(" << pk.x << "," << pk.y << "," << pk.z << ")";
    cout << " ftot=(" << ftot.x << "," << ftot.y << "," << ftot.z << ")" << endl;
}

// volume for contact calculations
void MatVelocityField::AddContactVolume(double vol) { volume += vol; }
void MatVelocityField::SetContactVolume(double vol) { volume = vol; }
double MatVelocityField::GetContactVolume(void) const { return volume; }

// Get velocity
Vector MatVelocityField::GetVelocity(void) { return vk; }

// moment zero for direction of velocity
// Let pk = (pk.norm) norm + (pk.tang) tang
// Here we want to subtract component in norm direction using pk - (pk.norm) norm
void MatVelocityField::SetMomentVelocityDirection(Vector *norm,int passType)
{	double dotn = DotVectors(&pk, norm);
	AddScaledVector(&pk, norm, -dotn);
	if(passType==UPDATE_MOMENTUM_CALL)
		AddScaledVector(&ftot,norm,-dotn/timestep);
}

// add moment for one component only
void MatVelocityField::AddMomentVelocityDirection(Vector *norm,double vel,int passType)
{	AddScaledVector(&pk, norm, mass*vel);
	if(passType==UPDATE_MOMENTUM_CALL)
		AddScaledVector(&ftot,norm,mass*vel/timestep);
}

// set symmetry plane momenta to zero
// only called when ADJUST_EXTRAPOLATED_PK_FOR_SYMMETRY is defined
void MatVelocityField::AdjustForSymmetryBC(int fixedDirection)
{   if(fixedDirection&XSYMMETRYPLANE_DIRECTION)
    {   pk.x = 0.;
		xpic[pkCopy].x=0.;
    }
    if(fixedDirection&YSYMMETRYPLANE_DIRECTION)
    {   pk.y = 0.;
		xpic[pkCopy].y=0.;
    }
    if(fixedDirection&ZSYMMETRYPLANE_DIRECTION)
    {   pk.z = 0.;
		xpic[pkCopy].z=0.;
    }
}

// Start adjusting force for nodal velocity BC
//      Ftot = (Ftot.norm) norm + (Ftot.tang) tang
// But, now we want it to start with
//      Ftotnew = -(pk.norm)/deltime norm + (Ftot.tang) tang
//      Ftotnew = Ftot - ((Ftot.norm) + (pk.norm)/deltime) norm
// Note that if have same norm on same node, the net result after first pass is
//      Ftot1 = -(pk.norm)/deltime norm + (Ftot.tang) tang
// Then on second pass, change in force is:
//		dotf = Ftot1.norm = -(pk.norm)/deltime
//      -dotf-dotp/deltime = -(pk.norm)/deltime + (pk.norm)/deltime = 0
// But might have physical issues if components of norm overlap on the same node such as
//    x axis and  skew xy or xz on the same node
void MatVelocityField::SetFtotDirection(Vector *norm,double deltime,Vector *freaction)
{   double dotf = DotVectors(&ftot, norm);
	double dotp = DotVectors(&pk, norm);
	// the change in force is (-dotf-dotp/deltime) norm, which is zero for second BC with same norm
	CopyScaleVector(freaction,norm,-dotf-dotp/deltime);
	AddVector(&ftot,freaction);
}

// add one component of force such that updated momentum will be mass*velocity
void MatVelocityField::AddFtotDirection(Vector *norm,double deltime,double vel,Vector *freaction)
{	// the change in force is (mass*vel/deltime) norm
	Vector deltaF;
	CopyScaleVector(&deltaF,norm,mass*vel/deltime);
	AddVector(&ftot,&deltaF);
	AddVector(freaction,&deltaF);
}

// total force
Vector MatVelocityField::GetFtot(void) const { return ftot; }
Vector *MatVelocityField::GetFtotPtr(void) { return &ftot; }

// rigid field
bool MatVelocityField::IsRigidField(void) const { return (bool)(flags&RIGID_FIELD_BIT); }
void MatVelocityField::SetRigidField(bool setting)
{	if(setting)
		flags |= RIGID_FIELD_BIT;
	else if(flags&RIGID_FIELD_BIT)
		flags -= RIGID_FIELD_BIT;
}

// Rigid contact, but not rigid block
bool MatVelocityField::IsFixedRigidField(void) const
{	if(!(flags&RIGID_FIELD_BIT)) return false;
	if(flags&RIGID_BLOCK_BIT) return false;
	return true;
}

// rigid block field
bool MatVelocityField::IsRigidBlockField(void) const { return (bool)(flags&RIGID_BLOCK_BIT); }

// all flags
int MatVelocityField::GetFlags() const { return flags; }

// ignoring cracks
bool MatVelocityField::IgnoresCracks(void) const
{	return (flags&IGORE_CRACKS_BIT) != 0 ? true : false;
}

#pragma mark CLASS METHODS

// return true if references field is active in this time step
bool MatVelocityField::ActiveField(MatVelocityField *mvf)
{ return mvf==NULL ? false : (mvf->numberPoints>0) ; }

// return true if references field that is active and is not for rigid materials
bool MatVelocityField::ActiveRigidField(MatVelocityField *mvf)
{ return mvf==NULL ? false : (mvf->numberPoints>0 && mvf->IsRigidField()) ; }

/* -------- Following look for non rigid fields of various types ----------- */

// return true if references field that is active and is not for rigid materials
bool MatVelocityField::ActiveNonrigidField(MatVelocityField *mvf)
{ return mvf==NULL ? false : (mvf->numberPoints>0 && !mvf->IsRigidField()) ; }

// return true if references nonrigid field in memory
//    if non-rigid material is ignoring cracks, then only true for field [0]
// In NairnMPM, all non rigid are source fields
bool MatVelocityField::ActiveNonrigidSourceField(MatVelocityField *mvf,int fieldNum)
{	return ActiveNonrigidField(mvf);
}

// return true if references field that is active, is not rigid, and matches requireCracks
// (i.e., if requireCracks is true, this field must see cracks, otherwise any field is OK)
// In NairnMPM, same as having any nonrigid particles
bool MatVelocityField::ActiveNonrigidSeesCracksField(MatVelocityField *mvf,bool requireCracks)
{	return ActiveNonrigidField(mvf);
}


