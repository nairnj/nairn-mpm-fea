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

// Add to force spread out over the materials so each has same extra accelerations = f/M_i
// Only called by AddTractionForce() and CrackInterfaceForce()
void CrackVelocityFieldSingle::AddFtotSpreadTask3(Vector *f)
{	if(MatVelocityField::ActiveField(mvf[0]))
		mvf[0]->AddFtot(f);
}

// Add grid dampiong force at a node in g mm/sec^2
void CrackVelocityFieldSingle::AddGridDampingTask3(double extDamping)
{	if(MatVelocityField::ActiveField(mvf[0]))
		mvf[0]->AddGridDampingTask3(extDamping);
}

#pragma mark TASK 4 METHODS

// update momenta for this MPM step
//  pk(i+1) = pk(i) + ftot * dt
void CrackVelocityFieldSingle::UpdateMomentaOnField(double timestep)
{	if(MatVelocityField::ActiveField(mvf[0]))
        mvf[0]->UpdateMomentum(timestep);
}

#pragma mark TASK 6 METHODS

// zero momentum and displacement at a node for new calculations
void CrackVelocityFieldSingle::RezeroNodeTask6(double)
{	if(MatVelocityField::ActiveField(mvf[0]))
	{	ZeroVector(&mvf[0]->pk);
		ZeroVector(&mvf[0]->disp);
		mvf[0]->SetContactVolume(0.);
	}
}

#pragma mark VELOCITY METHODS

// Calculate velocity at a node from current momentum and mass matrix in all velocity fields
void CrackVelocityFieldSingle::CalcVelocityForStrainUpdate(void)
{	if(MatVelocityField::ActiveField(mvf[0]))
		mvf[0]->CalcVelocityForStrainUpdate();
}

#pragma mark BOUNDARY CONDITIONS

// zero one component of moment and velocity
void CrackVelocityFieldSingle::SetMomVel(Vector *norm)
{	if(MatVelocityField::ActiveField(mvf[0]))
        mvf[0]->SetMomentVelocityDirection(norm);
}

// add one component of momentum and velocity from BCs
void CrackVelocityFieldSingle::AddMomVel(Vector *norm,double vel)
{	if(MatVelocityField::ActiveField(mvf[0]))
        mvf[0]->AddMomentVelocityDirection(norm,vel);
}

// set one component of force to -p(interpolated)/time such that updated momentum
//    of pk.i + deltime*ftot.i will be zero
void CrackVelocityFieldSingle::SetFtotDirection(Vector *norm,double deltime)
{	if(MatVelocityField::ActiveField(mvf[0]))
        mvf[0]->SetFtotDirection(norm,deltime);
}

// add one component of force such that updated momentum will be mass*velocity
void CrackVelocityFieldSingle::AddFtotDirection(Vector *norm,double deltime,double vel)
{	if(MatVelocityField::ActiveField(mvf[0]))
        mvf[0]->AddFtotDirection(norm,deltime,vel);
}

#pragma mark ACCESSORS

// total mass all velocity fields
double CrackVelocityFieldSingle::GetTotalMass(void)
{	return MatVelocityField::ActiveField(mvf[0]) ? mvf[0]->mass : 0. ;
}

// total mass and kinetric energy all velocity fields (rigid particles not counted)
// in g-mm2/sec^2 = nanoJ
void CrackVelocityFieldSingle::AddKineticEnergyAndMass(double &kineticEnergy,double &totalMass)
{	if(MatVelocityField::ActiveNonrigidField(mvf[0]))
	{	if(mvf[0]->mass > 0.)
		{	totalMass += mvf[0]->mass;
			double magp = DotVectors(&mvf[0]->pk,&mvf[0]->pk);
			kineticEnergy += 0.5*magp/mvf[0]->mass;
		}
	}
}

// get volume when only a single material (overridden when might be more)
double CrackVelocityFieldSingle::GetVolumeNonrigid(void)
{	return MatVelocityField::ActiveField(mvf[0]) ? mvf[0]->GetContactVolume() : 0. ;
}

// get volume when only a single material (overridden when might be more)
double CrackVelocityFieldSingle::GetVolumeTotal(double ndr)
{	return MatVelocityField::ActiveField(mvf[0]) ? mvf[0]->GetContactVolume() : 0. ;
}

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
        return mvf[0]->GetFtot();
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
	
