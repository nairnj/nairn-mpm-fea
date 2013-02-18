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
void MatVelocityField::IncrementNodalVelAcc(double fi,Vector *delv,Vector *dela)
{
    double mnode = fi/mass;					// Ni/mass
	AddScaledVector(delv,&pk,mnode);		// velocity += p/mass
	AddScaledVector(dela,&ftot,mnode);		// acceleration += f/mass
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

// moment and velocity zero for on component only
void MatVelocityField::SetMomentVelocityDirection(int dir)
{   if(dir==X_DIRECTION)
    {	pk.x = 0.;
        vk.x = 0.;
    }
    else if(dir==Y_DIRECTION)
    {	pk.y = 0.;
        vk.y = 0.;
    }
    else
    {	pk.z = 0.;
        vk.z = 0.;
    }
}

// add moment and velocity for one component only
void MatVelocityField::AddMomentVelocityDirection(int dir,double vel)
{   if(dir==X_DIRECTION)
    {	pk.x += mass*vel;
        vk.x += vel;
    }
    else if(dir==Y_DIRECTION)
    {	pk.y += mass*vel;
        vk.y += vel;
    }
    else
    {	pk.z += mass*vel;
        vk.z += vel;
    }
}

// set component of Ftot to -p/dt (used by boundary conditions)
void MatVelocityField::SetFtotDirection(int dir,double deltime)
{   if(dir==X_DIRECTION)
        ftot.x = -pk.x/deltime;
    else if(dir==Y_DIRECTION)
        ftot.y = -pk.y/deltime;
    else
        ftot.z = -pk.z/deltime;
}

// add one component of force such that updated momentum will be mass*velocity
void MatVelocityField::AddFtotDirection(int dir,double deltime,double vel)
{   if(dir==X_DIRECTION)
        ftot.x += mass*vel/deltime;
    else if(dir==Y_DIRECTION)
        ftot.y += mass*vel/deltime;
    else
        ftot.z += mass*vel/deltime;
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

