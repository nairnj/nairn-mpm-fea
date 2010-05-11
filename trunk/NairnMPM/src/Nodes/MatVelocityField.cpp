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

#pragma mark INITIALIZATION

// Constructors
MatVelocityField::MatVelocityField()
{	if(fmobj->multiMaterialMode)
		massGrad=new Vector;
	else
		massGrad=NULL;
	Zero();
}

// Destructor
MatVelocityField::~MatVelocityField()
{	if(massGrad!=NULL)
		delete massGrad;
}

// zero data at start of time step
void MatVelocityField::Zero(void)
{	mass=0.;
	ZeroVector(&pk);
	ZeroVector(&fint);
	ZeroVector(&fext);
	ZeroVector(&vk);
	ZeroVector(&ftot);
	ZeroVector(&disp);
	if(massGrad!=NULL) ZeroVector(massGrad);
	numberPoints=0;
}

#pragma mark METHODS

// add to momentum in task one. Save velocity for the first point in case only
// has one point (for more accurate velocity calculation)
void MatVelocityField::AddMomentumTask1(Vector *addPk,Vector *vel)
{	AddVector(&pk,addPk);
	numberPoints++;
	if(numberPoints==1) vk=*vel;
}

// in response to contact, change the momentum
// for a single point, calculate the velocity
// if postUpdate is true, then adjust ftot as well to keep in sync with momentum change
void MatVelocityField::ChangeMatMomentum(Vector *delP,bool postUpdate,double deltime)
{	AddVector(&pk,delP);
	if(numberPoints==1) CopyScaleVector(&vk,&pk,1./mass);
	if(postUpdate) AddScaledVector(&ftot, delP, 1./deltime);
}

// Calculate velocity at a node from current momentum and mass matrix in all velocity fields
void MatVelocityField::CalcVelocityForStrainUpdate(void)
{	// only 1 point or rigid contact material is stored already, 0 will have zero velocity
	if(numberPoints<=1 || mass==0) return;
	CopyScaleVector(&vk,&pk,1./mass);
}

// Calculate total force at a node from current values
// 		now m*a in g mm/sec^2
void MatVelocityField::CalcFtotTask3(double extDamping)
{	ftot.x=fint.x+fext.x-extDamping*pk.x;
	ftot.y=fint.y+fext.y-extDamping*pk.y;
	ftot.z=fint.z+fext.z-extDamping*pk.z;
}

#pragma mark ACCESSORS

// for debugging
void MatVelocityField::Describe(void)
{
	cout << "#      - Material Field: n="<<  numberPoints << " mass=" << mass << endl;
}

#pragma mark CLASS METHODS

// return true if references field is active in this time step
bool MatVelocityField::ActiveField(MatVelocityField *mvf) { return mvf==NULL ? (bool)FALSE : (mvf->numberPoints>0) ; }

// return true if references field that is active and is not for rigid materials
bool MatVelocityField::ActiveNonrigidField(MatVelocityField *mvf) { return mvf==NULL ? (bool)FALSE : (mvf->numberPoints>0 && mvf->mass>0.) ; }

// return true if references field that is active and is not for rigid materials
bool MatVelocityField::ActiveRigidField(MatVelocityField *mvf) { return mvf==NULL ? (bool)FALSE : (mvf->numberPoints>0 && mvf->mass<=0.) ; }

