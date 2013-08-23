/********************************************************************************
    MatVelocityField.cpp
    nairn-mpm-fea
    
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
#include "Cracks/CrackHeader.hpp"
#include "Nodes/NodalPoint.hpp"

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
	{	ZeroVector(&ftot);
#ifdef USE_FEXT
		ZeroVector(&fext);
#endif
	}
	ZeroVector(&vk);
	ZeroVector(&disp);
	if(volumeGrad!=NULL) ZeroVector(volumeGrad);
	numberPoints=0;
	volume=0.;
}

#pragma mark METHODS

// add to momentum in task 1. Save velocity for the first point in case only
// has one point (for more accurate velocity calculation)
// Momentum is g-mm/sec
void MatVelocityField::AddMomentumTask1(Vector *addPk,Vector *vel,int numPts)
{	AddVector(&pk,addPk);
	numberPoints += numPts;
	if(numberPoints==1) vk=*vel;
}

// copy mass and momentum from ghost to real node
void MatVelocityField::CopyMassAndMomentum(NodalPoint *real,int vfld,int matfld)
{
	// skip if none
	if(numberPoints==0) return;
	
	// momentum
	if(numberPoints==1)
		real->AddMomentumTask1((short)vfld,matfld,mass,&vk,1);
	else
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
	if(!rigidField)
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
	if(numberPoints==0 || rigidField) return;
	
	// momentum
	if(numberPoints==1)
		real->AddMomentumTask6((short)vfld,matfld,mass,&vk);
	else
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
	if(numberPoints==0) return;
	
	// Ftot
	real->AddFtotTask3(vfld,matfld,&ftot);
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
{	ftot.x += delP->x;
	ftot.y += delP->y;
	ftot.z += delP->z;
}

// Calculate velocity at a node from current momentum and mass matrix in all velocity fields
// Velocity is mm/sec
void MatVelocityField::CalcVelocityForStrainUpdate(void)
{	// only 1 point or rigid contact material is stored already, 0 will have zero velocity
	if(numberPoints<=1 || rigidField || mass==0.) return;
	CopyScaleVector(&vk,&pk,1./mass);
}

// Add grid dampiong force at a node in g mm/sec^2
void MatVelocityField::AddGridDampingTask3(double extDamping)
{	ftot.x -= extDamping*pk.x;
	ftot.y -= extDamping*pk.y;
	ftot.z -= extDamping*pk.z;
#ifdef USE_FEXT
	ftot.x += fext.x;
	ftot.y += fext.y;
	ftot.z += fext.z;
#endif
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
void MatVelocityField::Describe(int fldnum) const
{
	cout << "#      " << fldnum << ". Material Field: n="<<  numberPoints << " mass=" << mass ;
    cout << " p=(" << pk.x << "," << pk.y << "," << pk.z << ")" << endl;
}

// volume for contact calculations
void MatVelocityField::AddContactVolume(double vol) { volume += vol; }
void MatVelocityField::SetContactVolume(double vol) { volume = vol; }
double MatVelocityField::GetContactVolume(void) const { return volume; }

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
// Ftot = (Ftot.norm) norm + (Ftot.tang) tang, but now we want
// Ftotnew = -(pk.norm)/deltime norm + (Ftot.tang) tang
// Ftotnew = Ftot - (Ftot.norm) norm - (pk.norm)/deltime norm
// The reaction due to BC is -(pk.norm)/deltime norm
//
// Note that is have same norm on same node, the net result have first pass is
//    Ftot1 = -(pk.norm)/deltime norm + (Ftot.tang) tang
// Then on second pass, dotf = -(pk.norm)/deltime to give same result
//    Ftot2 = (-(pk.norm)/deltime+(pk.norm)/deltime-(pk.norm)/deltime) norm + (Ftot.tang) tang
//          = -(pk.norm)/deltime norm + (Ftot.tang) tang
// But might have physical issues if component of norm overlap on the save node such as
//    x axis ans  skew xy or xz on the same node
void MatVelocityField::SetFtotDirection(Vector *norm,double deltime,Vector *freaction)
{   double dotf = DotVectors(&ftot, norm);
	double dotp = DotVectors(&pk, norm);
	// the change in force is (-dotf - dotp/deltime) norm, which is zero for second BC with same norm
	CopyScaleVector(freaction,norm,-dotf - dotp/deltime);
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
Vector MatVelocityField::GetFtot(void) { return ftot; }
Vector *MatVelocityField::GetFtotPtr(void) { return &ftot; }

#pragma mark CLASS METHODS

// return true if references field is active in this time step
bool MatVelocityField::ActiveField(MatVelocityField *mvf)
{ return mvf==NULL ? (bool)FALSE : (mvf->numberPoints>0) ; }

// return true if references field that is active and is not for rigid materials
bool MatVelocityField::ActiveNonrigidField(MatVelocityField *mvf)
{ return mvf==NULL ? (bool)FALSE : (mvf->numberPoints>0 && !mvf->rigidField) ; }

// return true if references field that is active and is not for rigid materials
bool MatVelocityField::ActiveRigidField(MatVelocityField *mvf)
{ return mvf==NULL ? (bool)FALSE : (mvf->numberPoints>0 && mvf->rigidField) ; }

