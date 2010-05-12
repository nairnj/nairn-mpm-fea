/********************************************************************************
	CrackVelocityFieldSingle.cpp
	NairnMPM

	Created by John Nairn on 21 August 2009.
	Copyright (c) 2009 John A. Nairn, All rights reserved.
********************************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Nodes/CrackVelocityFieldSingle.hpp"
#include "Exceptions/CommonException.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"

#pragma mark INITIALIZATION

// constructor
CrackVelocityFieldSingle::CrackVelocityFieldSingle(short theLoc,int cnum) : CrackVelocityField(theLoc,cnum)
{
}

// zero the one material if it was active
void CrackVelocityFieldSingle::ZeroMatFields(void)
{	if(MatVelocityField::ActiveField(mvf[0]))
		mvf[0]->Zero();
}

#pragma mark TASK 1 METHODS

// Delete empty velocity fields, count number of materials, and return total mass
double CrackVelocityFieldSingle::GetTotalMassAndCount(void)
{	return MatVelocityField::ActiveField(mvf[0]) ? mvf[0]->mass : 0. ;
}

#pragma mark TASK 3 METHODS

// Add to fint spread out over the materials to each has same extra accelerations = f/M_i
// only called to add interface foces on a crack
void CrackVelocityFieldSingle::AddFintSpreadTask3(Vector *f)
{	if(MatVelocityField::ActiveField(mvf[0]))
		AddVector(&mvf[0]->fint,f);
}

// Add to fext spread out over the materials to each has same extra accelerations = f/M_i
// Only called for crack traction forces
void CrackVelocityFieldSingle::AddFextSpreadTask3(Vector *f)
{	if(MatVelocityField::ActiveField(mvf[0]))
		AddVector(&mvf[0]->fext,f);
}

// Calculate total force at a node from current values
// 		now m*a in g mm/sec^2
void CrackVelocityFieldSingle::CalcFtotTask3(double extDamping)
{	if(MatVelocityField::ActiveField(mvf[0]))
		mvf[0]->CalcFtotTask3(extDamping);
}

#pragma mark TASK 4 METHODS

// update momenta for this MPM step
//  pk(i+1) = pk(i) + ftot * dt
void CrackVelocityFieldSingle::UpdateMomentaTask4(double timestep)
{	if(MatVelocityField::ActiveField(mvf[0]))
		AddScaledVector(&mvf[0]->pk,&mvf[0]->ftot,timestep);
}

#pragma mark TASK 6 METHODS

// zero momentum and displacement at a node for new calculations
void CrackVelocityFieldSingle::RezeroNodeTask6(double)
{	if(MatVelocityField::ActiveField(mvf[0]))
	{	ZeroVector(&mvf[0]->pk);
		ZeroVector(&mvf[0]->disp);
	}
}

#pragma mark VELOCITY METHODS

// Calculate velocity at a node from current momentum and mass matrix in all velocity fields
void CrackVelocityFieldSingle::CalcVelocityForStrainUpdate(void)
{	if(MatVelocityField::ActiveField(mvf[0]))
		mvf[0]->CalcVelocityForStrainUpdate();
}

#pragma mark BOUNDARY CONDITIONS

// zero x moment and velocity
void CrackVelocityFieldSingle::SetXMomVel(void)
{	if(MatVelocityField::ActiveField(mvf[0]))
	{	mvf[0]->pk.x=0.;
		mvf[0]->vk.x=0.;
	}
}

// zero y moment and velocity
void CrackVelocityFieldSingle::SetYMomVel(void)
{	if(MatVelocityField::ActiveField(mvf[0]))
	{	mvf[0]->pk.y=0.;
		mvf[0]->vk.y=0.;
	}
}

// zero z moment and velocity
void CrackVelocityFieldSingle::SetZMomVel(void)
{	if(MatVelocityField::ActiveField(mvf[0]))
	{	mvf[0]->pk.z=0.;
		mvf[0]->vk.z=0.;
	}
}

// zero momentum in direction (cos(angle), -sin(angle)) or vector rotated from postive x axis
// by clockwise angle. The desired vector is (p.t)t where t is unit vector normal to
// skew direction and here t = (sin(angle), cos(angle))
void CrackVelocityFieldSingle::SetSkewMomVel(double angle)
{	if(MatVelocityField::ActiveField(mvf[0]))
	{	double c=cos(angle),s=sin(angle);
		Vector *npk=&mvf[0]->pk;
		double momx=npk->x*s*s + npk->y*c*s;
		double momy=npk->x*c*s + npk->y*c*c;
		npk->x=momx;
		npk->y=momy;
		mvf[0]->vk.x=momx/mvf[0]->mass;
		mvf[0]->vk.y=momy/mvf[0]->mass;
	}
}

// add to x momentum and velocity from BCs
void CrackVelocityFieldSingle::AddXMomVel(double vx)
{	if(MatVelocityField::ActiveField(mvf[0]))
	{	mvf[0]->pk.x+=mvf[0]->mass*vx;
		mvf[0]->vk.x+=vx;
	}
}

// add to y momentum and velocity from BCs
void CrackVelocityFieldSingle::AddYMomVel(double vy)
{	if(MatVelocityField::ActiveField(mvf[0]))
	{	mvf[0]->pk.y+=mvf[0]->mass*vy;
		mvf[0]->vk.y+=vy;
	}
}

// add to z momentum and velocity from BCs
void CrackVelocityFieldSingle::AddZMomVel(double vz)
{	if(MatVelocityField::ActiveField(mvf[0]))
	{	mvf[0]->pk.z+=mvf[0]->mass*vz;
		mvf[0]->vk.z+=vz;
	}
}

// Add velocity in direction (cos(angle), -sin(angle)) or vector rotated from postive x axis
// by clockwise angle.
void CrackVelocityFieldSingle::AddSkewMomVel(double vel,double angle)
{	if(MatVelocityField::ActiveField(mvf[0]))
	{	double velx=cos(angle)*vel;
		double vely=-sin(angle)*vel;
		mvf[0]->pk.x+=mvf[0]->mass*velx;
		mvf[0]->vk.x+=velx;
		mvf[0]->pk.y+=mvf[0]->mass*vely;
		mvf[0]->vk.y+=vely;
	}
}

// set x force to -p(interpolated)/time such that updated momentum
//    of pk.x + deltime*ftot.x will be zero
void CrackVelocityFieldSingle::SetXFtot(double deltime)
{	if(MatVelocityField::ActiveField(mvf[0]))
		mvf[0]->ftot.x=-mvf[0]->pk.x/deltime;
}

// set y force to -p(interpolated)/time such that updated momentum
//    of pk.y + deltime*ftot.y will be zero
void CrackVelocityFieldSingle::SetYFtot(double deltime)
{	if(MatVelocityField::ActiveField(mvf[0]))
		mvf[0]->ftot.y=-mvf[0]->pk.y/deltime;
}

// set z force to -p(interpolated)/time such that updated momentum
//    of pk.z + deltime*ftot.z will be zero
void CrackVelocityFieldSingle::SetZFtot(double deltime)
{	if(MatVelocityField::ActiveField(mvf[0]))
		mvf[0]->ftot.z=-mvf[0]->pk.z/deltime;
}

// Change current ftot such that updated momentum of (pk.x + deltime*ftot.x, pk.y + deltime*ftot.y) will have zero
// momentum in the (cos(theta), -sin(angle)) direction (or direction clockwise from positive x axis by angle).
// Superpose force to induce zero momentum in the skew direction:
//			f dt = (p.t)t - p(interpolated)
// where t = (sin(angle),cos(angle)) is tangential vector, with the existing component of total force
// in the t direction or (f.t)t
void CrackVelocityFieldSingle::SetSkewFtot(double deltime,double angle)
{   if(MatVelocityField::ActiveField(mvf[0]))
	{	double c=cos(angle),s=sin(angle);
		double dfxdt=-mvf[0]->pk.x*c*c + mvf[0]->pk.y*c*s;	// to get zero skew momentum
		double dfydt=mvf[0]->pk.x*c*s - mvf[0]->pk.y*s*s;
		double fx=mvf[0]->ftot.x*s*s + mvf[0]->ftot.y*c*s;	// f normal to skew direction
		double fy=mvf[0]->ftot.x*c*s + mvf[0]->ftot.y*c*c;
		mvf[0]->ftot.x=fx + dfxdt/deltime;
		mvf[0]->ftot.y=fy + dfydt/deltime;
	}
}

// add to x force such that updated momentum will be mass*velocity
void CrackVelocityFieldSingle::AddXFtot(double deltime,double velx)
{	if(MatVelocityField::ActiveField(mvf[0]))
		mvf[0]->ftot.x+=mvf[0]->mass*velx/deltime;
}

// add to y force such that updated momentum will be mass*velocity
void CrackVelocityFieldSingle::AddYFtot(double deltime,double vely)
{	if(MatVelocityField::ActiveField(mvf[0]))
		mvf[0]->ftot.y+=mvf[0]->mass*vely/deltime;
}

// add to z force such that updated momentum will be mass*velocity
void CrackVelocityFieldSingle::AddZFtot(double deltime,double velz)
{	if(MatVelocityField::ActiveField(mvf[0]))
		mvf[0]->ftot.z+=mvf[0]->mass*velz/deltime;
}

// set force in the skew direction (cos(angle),-sin(angle)), or direction rotated from postive
// x axis by clockwise angle.
void CrackVelocityFieldSingle::AddSkewFtot(double deltime,double vel,double angle)
{	if(MatVelocityField::ActiveField(mvf[0]))
	{	double velx=cos(angle)*vel;
		double vely=-sin(angle)*vel;
		mvf[0]->ftot.x+=mvf[0]->mass*velx/deltime;
		mvf[0]->ftot.y+=mvf[0]->mass*vely/deltime;
	}
}

#pragma mark ACCESSORS

// total mass all velocity fields
double CrackVelocityFieldSingle::GetTotalMass(void)
{	return MatVelocityField::ActiveField(mvf[0]) ? mvf[0]->mass : 0. ;
}

// total mass all velocity fields
double CrackVelocityFieldSingle::GetMass(int matfld) { return GetTotalMass(); }

// get center of mass momentum for all material fields in this crack velocity field
Vector CrackVelocityFieldSingle::GetCMatMomentum(void)
{	if(MatVelocityField::ActiveField(mvf[0]))
		return mvf[0]->pk;
	else
	{	Vector pk;
		ZeroVector(&pk);
		return pk;
	}
}

// get center of mass displacement (actually sum of displacement*mass so displacement is vector/total mass)
Vector CrackVelocityFieldSingle::GetCMDisplacement(void)
{	if(MatVelocityField::ActiveField(mvf[0]))
		return mvf[0]->disp;
	else
	{	Vector dk;
		ZeroVector(&dk);
		return dk;
	}
}

// get center of mass momentum for all material fields in this crack velocity field
Vector CrackVelocityFieldSingle::GetCMatFtot(void)
{	if(MatVelocityField::ActiveField(mvf[0]))
		return mvf[0]->ftot;
	else
	{	Vector fk;
		ZeroVector(&fk);
		return fk;
	}
}

/* in response to crack contact, change moment by changing velocity of all materials
	the same amount
 
   Change velocity by dP/M, where M is total mass
   Material i velocity becomes vi = pi/mi + dP/M
   Material i momentum change is mi vi = pi + mi dP/M
*/
void CrackVelocityFieldSingle::ChangeMomentum(Vector *delP,bool postUpdate,double deltime)
{	if(MatVelocityField::ActiveField(mvf[0]))
		mvf[0]->ChangeMatMomentum(delP,postUpdate,deltime);
}

// copy all material velocity fields for boundary conditions methods, returning new offset into the save array
int CrackVelocityFieldSingle::CopyFieldMomenta(Vector *holdPk,int offset)
{	if(MatVelocityField::ActiveField(mvf[0]))
	{	holdPk[offset]=mvf[0]->pk;
		offset++;
	}
	return offset;
}

// paste all material velocity fields back for boundary conditions methods, returning new offset into the saved array
int CrackVelocityFieldSingle::PasteFieldMomenta(Vector *holdPk,int offset)
{	if(MatVelocityField::ActiveField(mvf[0]))
	{	mvf[0]->pk=holdPk[offset];
		offset++;
	}
	return offset;
}

// for debugging
void CrackVelocityFieldSingle::Describe(void)
{	CrackVelocityField::Describe();
	cout << "#     single material" << endl;
	if(MatVelocityField::ActiveField(mvf[0])) mvf[0]->Describe();
}
	
