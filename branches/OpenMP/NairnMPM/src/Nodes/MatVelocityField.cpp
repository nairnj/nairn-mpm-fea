/********************************************************************************
    MatVelocityField.cpp
    NairnMPM
    
    Created by John Nairn on 3 April 2009.
    Copyright (c) 2009 John A. Nairn, All rights reserved.
 
	Special case for rigid material
		pk will be sum of Vp*vel
		vk will be the rigid material current velocity
		disp will be sum Vp*disp
		mass grad will be volume gradient
		unscaled volume if the volume
	
********************************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/MatVelocityField.hpp"
#include "Exceptions/CommonException.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"

#pragma mark INITIALIZATION

// Constructors
MatVelocityField::MatVelocityField(short forRigid)
{	if(fmobj->multiMaterialMode)
		volumeGrad=new Vector;
	else
		volumeGrad=NULL;
	rigidField=FALSE;
	Zero();
	rigidField=forRigid;
}

// Destructor
MatVelocityField::~MatVelocityField()
{	if(volumeGrad!=NULL)
		delete volumeGrad;
}

// zero data at start of time step
void MatVelocityField::Zero(void)
{	mass=0.;
	ZeroVector(&pk);
	if(!rigidField)
	{	ZeroVector(&fint);
		ZeroVector(&fext);
		ZeroVector(&ftot);
	}
	ZeroVector(&vk);
	ZeroVector(&disp);
	if(volumeGrad!=NULL) ZeroVector(volumeGrad);
	numberPoints=0;
	volume=0.;
}

#pragma mark METHODS

// add to momentum in task one. Save velocity for the first point in case only
// has one point (for more accurate velocity calculation)
// Momentum is g-mm/sec
void MatVelocityField::AddMomentumTask1(Vector *addPk,Vector *vel)
{	AddVector(&pk,addPk);
	numberPoints++;
	if(numberPoints==1) vk=*vel;
}

// in response to contact, change the momentum
// for a single point, calculate the velocity
// if postUpdate is true, then adjust ftot as well to keep in sync with momentum change
// Momentum is g-mm/sec
void MatVelocityField::ChangeMatMomentum(Vector *delP,bool postUpdate,double deltime)
{	AddVector(&pk,delP);
	if(numberPoints==1) CopyScaleVector(&vk,&pk,1./mass);
	if(postUpdate) AddScaledVector(&ftot, delP, 1./deltime);
}

// in response to contact, and only for rigid materials, add contact force to ftot
// Contact force is g-mm/sec^2 = micro N (because N is kg-m/sec^2)
// Here just sum momenta applied to non-rigid materials. To get forces, this result
//	is multiplied by -1/timestep
void MatVelocityField::AddContactForce(Vector *delP)
{	AddVector(&ftot,delP);
}

// Calculate velocity at a node from current momentum and mass matrix in all velocity fields
// Velocity is mm/sec
void MatVelocityField::CalcVelocityForStrainUpdate(void)
{	// only 1 point or rigid contact material is stored already, 0 will have zero velocity
	if(numberPoints<=1 || rigidField || mass==0.) return;
	CopyScaleVector(&vk,&pk,1./mass);
}

// Calculate total force at a node from current values
// Now force = m*a in g mm/sec^2 = micro N (because N is kg-m/sec^2)
void MatVelocityField::CalcFtotTask3(double extDamping)
{	ftot.x = fint.x + fext.x - extDamping*pk.x;
	ftot.y = fint.y + fext.y - extDamping*pk.y;
	ftot.z = fint.z + fext.z - extDamping*pk.z;
}

// internal force - add or scale and add
void MatVelocityField::AddFint(Vector *f) { AddVector(&fint,f); }
void MatVelocityField::AddFint(Vector *f,double scaled) { AddScaledVector(&fint,f,scaled); }

// internal force - add or scale and add
void MatVelocityField::AddFext(Vector *f) { AddVector(&fext,f); }
void MatVelocityField::AddFext(Vector *f,double scaled) { AddScaledVector(&fext,f,scaled); }

// Update momentum for this MPM step
//  pk(i+1) = pk(i) + ftot * dt
void MatVelocityField::UpdateMomentum(double timestep)
{	// update momenta
    AddScaledVector(&pk,&ftot,timestep);
}

// on strain updates, increment nodal velocity and acceleration
// fi is shape function
void MatVelocityField::IncrementNodalVelAcc(double fi,Vector *delv,Vector *dela) const
{
    double mnode = fi/mass;					// Ni/mass
	delv->x += pk.x*mnode;					// velocity += p/mass
	delv->y += pk.y*mnode;
	delv->z += pk.z*mnode;
	dela->x += ftot.x*mnode;				// acceleration += f/mass
	dela->y += ftot.y*mnode;
	dela->z += ftot.z*mnode;
}

#pragma mark ACCESSORS

// for debugging
void MatVelocityField::Describe(void)
{
	cout << "#      - Material Field: n="<<  numberPoints << " mass=" << mass << endl;
}

// volume for contact calculations
void MatVelocityField::AddContactVolume(double vol) { volume += vol; }
void MatVelocityField::SetContactVolume(double vol) { volume = vol; }
double MatVelocityField::GetContactVolume(void) { return volume; }

// velocity
void MatVelocityField::SetVelocity(Vector *vel) { vk = *vel; }
Vector MatVelocityField::GetVelocity(void) { return vk; }
Vector *MatVelocityField::GetVelocityPtr(void) { return &vk; }

// moment and velocity zero for direction of velocity
// Let p = (pn.norm) norm + (pt.tang) tang (or same for v)
// Here we want to remove pn.norm using p - (pn.norm) norm
void MatVelocityField::SetMomentVelocityDirection(Vector *norm)
{	double dotn = DotVectors(&pk, norm);
	AddScaledVector(&pk, norm, -dotn);
	dotn = DotVectors(&vk, norm);
	AddScaledVector(&vk, norm, -dotn);
}

// add moment and velocity for one component only
void MatVelocityField::AddMomentVelocityDirection(Vector *norm,double vel)
{	AddScaledVector(&pk, norm, mass*vel);
	AddScaledVector(&vk, norm, vel);
}

// Set component of Ftot to -p/dt in velocity BC direction (used by boundary conditions)
// Ftot = (Ftot.norm) norm + (Ftot.tang) tang, but not we want
// Ftotnew = -(pk.norm)/deltime norm + (Ftot.tang) tang
// Ftotnew = Ftot - (Ftot.norm) norm -(pk.norm)/deltime norm
void MatVelocityField::SetFtotDirection(Vector *norm,double deltime)
{   double dotf = DotVectors(&ftot, norm);
	double dotp = DotVectors(&pk, norm);
	AddScaledVector(&ftot, norm, -dotf - dotp/deltime);
}

// add one component of force such that updated momentum will be mass*velocity
void MatVelocityField::AddFtotDirection(Vector *norm,double deltime,double vel)
{   AddScaledVector(&ftot, norm, mass*vel/deltime);
}

// total force
Vector MatVelocityField::GetFtot(void) { return ftot; }
Vector *MatVelocityField::GetFtotPtr(void) { return &ftot; }

#pragma mark CLASS METHODS

// return true if references field is active in this time step
bool MatVelocityField::ActiveField(MatVelocityField *mvf) { return mvf==NULL ? (bool)FALSE : (mvf->numberPoints>0) ; }

// return true if references field that is active and is not for rigid materials
bool MatVelocityField::ActiveNonrigidField(MatVelocityField *mvf) { return mvf==NULL ? (bool)FALSE : (mvf->numberPoints>0 && !mvf->rigidField) ; }

// return true if references field that is active and is not for rigid materials
bool MatVelocityField::ActiveRigidField(MatVelocityField *mvf) { return mvf==NULL ? (bool)FALSE : (mvf->numberPoints>0 && mvf->rigidField) ; }

