/********************************************************************************
    MatVelocityField.cpp
    nairn-mpm-fea
    
    Created by John Nairn on 3 April 2009.
    Copyright (c) 2009 John A. Nairn, All rights reserved.
 
	Special case for rigid material
		pk will be sum of Vp*vel
 		ContactInfo
			disp will be sum Vp*(pos-origpos), expos will be sum of Vp*pos
			and volume will be unscaled volume
********************************************************************************/
#if defined ( _MSC_VER) || defined (__APPLE__) 
#include "stdafx.h"
#endif

#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/MatVelocityField.hpp"
#include "Exceptions/CommonException.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "System/MPMPrefix.hpp"


// class statics
int MatVelocityField::pkCopy=0;

#pragma mark INITIALIZATION

// Constructors
// throws std::bad_alloc
MatVelocityField::MatVelocityField(int setFlags)
{	// for contact
	if(fmobj->multiMaterialMode || firstCrack!=NULL)
	{	contactInfo = new ContactTerms;
		if(mpmgrid.numContactVectors>0)
			contactInfo->terms = new Vector[mpmgrid.numContactVectors];
		else
			contactInfo->terms = NULL;
	}
	else
		contactInfo = NULL;

	// XPIC in single material mode needs 3 vectors (v* and two working copies in XPIC tasks)
	//		FMPM needs just 2 (v-v*) in [0] and two working copies in XPIC tasks)
	// XPIC in multimaterial modes needs 6 vectors (v* four working copies for XPIC tasks, stored delta p^alpha) (for MM_XPIC)
	//		FMPM also needs 6, extra is for persistent store of delta p due to contact
	// ...note that DELTA_VSTAR_PREV must be VSTAR_PREV+2 - XPIC calcs assume that
	// FLIP and PIC need only 1 vector (to copy momenta)
	
	// add one more to store pk outside the xpic space
	int numVecs = bodyFrc.XPICVectors()+1;
	vk = new Vector[numVecs];
	pkCopy = numVecs-1;
	
	SetRigidField(false);
	Zero();
	flags = setFlags;
}

// Destructor
MatVelocityField::~MatVelocityField()
{	// delete contact data
	if(contactInfo!=NULL)
	{	if(contactInfo->terms!=NULL)
			delete contactInfo->terms;
		delete contactInfo;
	}
	delete [] vk;
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
	ZeroVector(vk);
	ZeroVector(&vk[pkCopy]);
	ZeroContactTerms();
	numberPoints=0;
}

#pragma mark METHODS

// add to momentum in task 1. Legacy momentum is g-mm/sec
void MatVelocityField::AddMomentumTask1(Vector *addPk,Vector *vel,int numPts)
{	AddPk(addPk);
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
	{	real->AddContactTerms((short)vfld,matfld,contactInfo);
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
	{	real->AddContactTerms((short)vfld,matfld,contactInfo);
	}
}

// Make copy of momenta, to be restored after strain update and forces
// Zero MM contact if needed
// return mass in this field
double MatVelocityField::GetTotalMassAndCount(void)
{
	// copy the extrapolated momenta
	vk[pkCopy] = pk;
	
#if MM_XPIC == 1
	// When doing XPIC in multimaterial modes, need to track change in momenta too
	if(bodyFrc.UsingVstar()==VSTAR_WITH_CONTACT)
		ZeroVector(&vk[DELTA_VSTORE_VEC]);
#endif
	
	return mass;
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
	pk = vk[pkCopy];
	
#ifdef CHECK_NAN
	if(pk.x!=pk.x || pk.y!=pk.y || pk.z!=pk.z || ftot.x!=ftot.x || ftot.y!=ftot.y || ftot.z!=ftot.z)
	{
#pragma omp critical (output)
		{	cout << "\n# MatVelocityField::RestoreMomenta: stored momenta or intial force was corrupted" << endl;
			cout << "#      It was stored in vk[" << pkCopy << "]" << endl;
			Describe(-1);
		}
	}
#endif
}

// in response to contact, change the momentum
// for a single point, calculate the velocity
// callType == MASS_MOMENTUM_CALL, UPDATE_MOMENTUM_CALL, UPDATE_STRAINS_LAST_CALL
void MatVelocityField::ChangeMatMomentum(Vector *delP,int callType,double deltime)
{
	// add to momentum
	AddPk(delP);
	
	// callType dependent changes
	if(callType==UPDATE_MOMENTUM_CALL)
	{	// add to force too
		AddFtotScaled(delP, 1./deltime);
	}
	else if(callType==MASS_MOMENTUM_CALL)
	{	// hack to skip the first contact corrections - i.e.,do not change anything
		//return;
		
		// for contact to be correct, need contact correction to be in initial momenta
		vk[pkCopy] = pk;
#if MM_XPIC == 1
		// For XPIC in multimaterial mode, track sum of contact changes
		if(bodyFrc.UsingVstar()==VSTAR_WITH_CONTACT)
			AddVector(&vk[DELTA_VSTORE_VEC],delP);
	}
	else
	{	// update strains last, for FMPM in multimaterial mode, track sum of contact changes in all modes
		if(bodyFrc.UsingVstar()==VSTAR_WITH_CONTACT)
		{	AddVector(&vk[DELTA_VSTORE_VEC],delP);
		}
	}
#else
	}
#endif
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

// Several calculations for grid values
// 1. Velocity prior to strain update
// Note: do not call for rigid fields in multimaterial mode
void MatVelocityField::GridValueCalculation(int calcOption)
{	// exit if nothing there
	if(numberPoints==0 || mass==0.) return;
	
	switch(calcOption)
	{	case VELOCITY_FOR_STRAIN_UPDATE:
			// get velocity
			CopyScaleVector(vk,&pk,1./mass);
			break;
		default:
			break;
	}
}

// Add grid damping force at a node in g mm/sec^2
void MatVelocityField::AddGravityAndBodyForceTask3(Vector *gridBodyForce)
{	ftot.x += mass*gridBodyForce->x;
	ftot.y += mass*gridBodyForce->y;
	ftot.z += mass*gridBodyForce->z;
}

#ifdef RESTART_OPTION
// chack if velocity and accelerations are too high
// distant traveled = (p + 0.5*f*dt)*dt/m = v*dt + 0.5*a*dt^2
bool MatVelocityField::IsTravelTooMuch(double dt,double maxDist) const
{
    Vector dx = pk;
    AddScaledVector(&dx,&ftot,0.5*dt);
    ScaleVector(&dx,dt/mass);
    double dist = sqrt(DotVectors(&dx,&dx));
    return dist>maxDist;
}
#endif

// total momentum - add
void MatVelocityField::AddPk(Vector *f)
{	pk.x += f->x;
	pk.y += f->y;
	pk.z += f->z;
}
// total momentum - scale and add
void MatVelocityField::AddPkScaled(Vector *f,double scaled)
{	pk.x += f->x*scaled;
	pk.y += f->y*scaled;
	pk.z += f->z*scaled;
}
// total force - add
void MatVelocityField::AddFtot(Vector *f)
{	ftot.x += f->x;
	ftot.y += f->y;
	ftot.z += f->z;
}
// internal force - scale and add
void MatVelocityField::AddFtotScaled(Vector *f,double scaled)
{	ftot.x += f->x*scaled;
	ftot.y += f->y*scaled;
	ftot.z += f->z*scaled;
}

// Update momentum for this MPM step
//  pk(i+1) = pk(i) + ftot * dt
void MatVelocityField::UpdateMomentum(double timestep)
{	// update momenta
    AddPkScaled(&ftot,timestep);
	
#ifdef CHECK_NAN
	if(pk.x!=pk.x || pk.y!=pk.y || pk.z!=pk.z)
	{
#pragma omp critical (output)
		{	cout << "\n# MatVelocityField::UpdateMomentum: updated momentum corrupted" << endl;
			Describe(-1);
		}
	}
#endif
}

// Support XPIC calculations
void MatVelocityField::XPICSupport(int xpicCalculation,int xpicOption,NodalPoint *real,double timestep,int m,int k,double vsign)
{
	switch(xpicCalculation)
	{	case INITIALIZE_XPIC:
		{	// these are needed zero on ghost and real nodes
			ZeroVector(&vk[VSTARNEXT_VEC]);
#if MM_XPIC == 1
			if(bodyFrc.UsingVstar()==VSTAR_WITH_CONTACT)
				ZeroVector(&vk[DELTA_VSTARNEXT_VEC]);
#endif
			
			// skip if none or if ghost node setting
			if(numberPoints==0 || timestep<0.) return;
			
			// set vStarPrev to k*v^+ = k*pi^+/mi, which is updated velocity
			double korder = (double)bodyFrc.GetXPICOrder();
			if(bodyFrc.UsingFMPM())
			{	CopyScaleVector(&vk[VSTARPREV_VEC],&pk,korder/mass);
				
				// set v* = vStarPrev (build velocity directly in v^*)
				vk[VSTAR_VEC] = vk[VSTARPREV_VEC];
			}
			else
			{	// For XPIC (only called during particle update), set to k*v = k(pi^+-fi*dt)/mi
				vk[VSTARPREV_VEC] = pk;
				AddScaledVector(&vk[VSTARPREV_VEC], &ftot, -timestep);
				ScaleVector(&vk[VSTARPREV_VEC],korder/mass);
				
				// set v* = vStarPrev * a*dt (build v-v^* + a*dt directly in v*)
				vk[VSTAR_VEC] = vk[VSTARPREV_VEC];
				AddScaledVector(&vk[VSTAR_VEC], &ftot, timestep/mass);
			}
			
			
#if MM_XPIC == 1
			// in multimaterial mode, get delta v from stored delta p due to contact
			if(bodyFrc.UsingVstar()==VSTAR_WITH_CONTACT)
			{	CopyScaleVector(&vk[DELTA_VSTARPREV_VEC],&vk[DELTA_VSTORE_VEC],(korder-1.)/mass);
			}
#endif
			break;
		}
		
		case UPDATE_VSTAR:
		{	// Add to Vstar and reset for next iteration
			if(numberPoints==0) return;
			
			// add to vStar += (-1)^k * vStarNext
			vk[VSTAR_VEC].x += vsign*vk[VSTARNEXT_VEC].x;
			vk[VSTAR_VEC].y += vsign*vk[VSTARNEXT_VEC].y;
			vk[VSTAR_VEC].z += vsign*vk[VSTARNEXT_VEC].z;
			
			// copy to previous
			vk[VSTARPREV_VEC] = vk[VSTARNEXT_VEC];
			
			// zero vStarNext
			ZeroVector(&vk[VSTARNEXT_VEC]);
			
#if MM_XPIC == 1
			if(bodyFrc.UsingVstar()==VSTAR_WITH_CONTACT)
			{	// add to v += (-1)^k * deltaVStarPrev
				// This adding previous
				vk[VSTAR_VEC].x += vsign*vk[DELTA_VSTARPREV_VEC].x;
				vk[VSTAR_VEC].y += vsign*vk[DELTA_VSTARPREV_VEC].y;
				vk[VSTAR_VEC].z += vsign*vk[DELTA_VSTARPREV_VEC].z;
				
				// copy to previous
				vk[DELTA_VSTARPREV_VEC] = vk[DELTA_VSTARNEXT_VEC];
				
				// zero deltaVStarNext
				ZeroVector(&vk[DELTA_VSTARNEXT_VEC]);
			}
#endif
			break;
		}
			
		case COPY_VSTARNEXT:
		{	// Copy to real node and zero vStarNext on this ghose node
			if(numberPoints==0) return;
			
			// add to real node (note that vfld in xpicOption and matfld in m
			Vector *vStarNextj = &vk[VSTARNEXT_VEC];
			real->AddVStarNext((short)xpicOption,m,vStarNextj,NULL,NULL,NULL,1.,1.);
			
			// zero on this material field on this ghost node for next loop
			ZeroVector(vStarNextj);
#if MM_XPIC == 1
			// Assume vector for multimaterial mode is shifted by 2
			if(bodyFrc.UsingVstar()==VSTAR_WITH_CONTACT) ZeroVector(vStarNextj+2);
#endif
			break;
		}
	
		default:
			break;
	}
}

// add to vStarNext
void MatVelocityField::AddVStarNext(Vector *vStarPrevj,Vector *delXiMpPtr,Vector *delXjPtr,
									  Matrix3 *Dpinv,double weight,double weightContact)
{
	// no particle spin
	vk[VSTARNEXT_VEC].x += weight*vStarPrevj->x;
	vk[VSTARNEXT_VEC].y += weight*vStarPrevj->y;
	vk[VSTARNEXT_VEC].z += weight*vStarPrevj->z;
	
#if MM_XPIC == 1
	// add contact term - assume in vector +2 from vStarPrevj
	if(bodyFrc.UsingVstar()==VSTAR_WITH_CONTACT)
	{	Vector *deltaVStarPrevj = vStarPrevj+2;
		vk[DELTA_VSTARNEXT_VEC].x += weightContact*deltaVStarPrevj->x;
		vk[DELTA_VSTARNEXT_VEC].y += weightContact*deltaVStarPrevj->y;
		vk[DELTA_VSTARNEXT_VEC].z += weightContact*deltaVStarPrevj->z;
	}
#endif
}

// on particle updates, increment nodal velocity and acceleration and others as needed
// fi is shape function
void MatVelocityField::IncrementNodalVelAcc(double fi,GridToParticleExtrap *gp) const
{
 	// Summing S v+ (FMPM has Sv+(k) and XPIC has Sv(k))
	gp->Svtilde.x += vk->x*fi;				// velocity += v[0]
	gp->Svtilde.y += vk->y*fi;
	gp->Svtilde.z += vk->z*fi;
	
	if(gp->m<=0)
	{	double mnode = fi/mass;					// Ni/mass
	
		// Summing S a (for FLIP and XPIC)
		gp->Sacc.x += ftot.x*mnode;
		gp->Sacc.y += ftot.y*mnode;
		gp->Sacc.z += ftot.z*mnode;
		
		if(gp->m<-1)
		{	// Summing S v^+ lumped for XPIC(k>1)
			gp->Svlumped.x += pk.x*mnode;
			gp->Svlumped.y += pk.y*mnode;
			gp->Svlumped.z += pk.z*mnode;
		}
	}
}

// zero momentum and displacement at a node for new calculations
void MatVelocityField::RezeroNodeTask6(void)
{	ZeroVector(&pk);
	ZeroContactTerms();
#if MM_XPIC == 1
	// When doing FMPM in multimaterial modes, rzero vector to track delta p
	if(bodyFrc.UsingVstar()==VSTAR_WITH_CONTACT)
		ZeroVector(&vk[DELTA_VSTORE_VEC]);
#endif
}

#pragma mark CONTACTINFO METHODS

// Zero contact terms (if used)
void MatVelocityField::ZeroContactTerms(void)
{	if(contactInfo==NULL) return;
	contactInfo->cvolume=0.;
	if(contactInfo->terms!=NULL)
	{	for(int i=0;i<mpmgrid.numContactVectors;i++)
			ZeroVector(&contactInfo->terms[i]);
	}
}

// pointer to extrapolated position or displacement
// must have contactInfo and must have the requested term
const Vector *MatVelocityField::GetContactDispPtr(bool useDisps) const
{	return useDisps ? &contactInfo->terms[mpmgrid.displacementIndex] :
						&contactInfo->terms[mpmgrid.positionIndex] ;
}

// add to a contact extrapolation vector with scaling. Do not call unless the index>=0
// must have contactInfo and must have the specified index
void MatVelocityField::AddContactVector(int index,Vector *toadd,double weight)
{	AddScaledVector(&contactInfo->terms[index],toadd,weight);
}

// add toa contact extrapolation vector. Do not call unless the index>=0
// must have contactInfo  and must have the specified index
void MatVelocityField::AddContactVector(int index,Vector *toadd)
{	AddVector(&contactInfo->terms[index],toadd);
}

// volume for contact calculations
void MatVelocityField::AddContactVolume(double vol) { contactInfo->cvolume += vol; }
void MatVelocityField::SetContactVolume(double vol) { contactInfo->cvolume = vol; }
double MatVelocityField::GetContactVolume(void) const { return contactInfo->cvolume; }

#pragma mark ACCESSORS

// for debugging
void MatVelocityField::Describe(int fldnum) const
{
	cout << "#      " << fldnum << ". Material Field: n="<<  numberPoints << " mass=" << mass ;
    cout << " p=(" << pk.x << "," << pk.y << "," << pk.z << ")";
    cout << " ftot=(" << ftot.x << "," << ftot.y << "," << ftot.z << ")" << endl;
}

// get velocity
Vector MatVelocityField::GetVelocity(void) { return vk[0]; }

// Get vStarPrev pointer
Vector *MatVelocityField::GetVStarPrev(void) const { return &vk[VSTARPREV_VEC]; };

// moment zero for direction of velocity
// Let pk = (pk.norm) norm + (pk.tang) tang
// Here we want to subtract component in norm direction using pk - (pk.norm) norm
// When have forces, subtract their normal component too
void MatVelocityField::ZeroVelocityBC(Vector *norm,int passType,double deltime,Vector *freaction)
{
	double dotn;
	
	switch(passType)
	{	case UPDATE_MOMENTUM_CALL:
			dotn = DotVectors(&pk, norm);
			AddPkScaled(norm, -dotn);
			AddFtotScaled(norm,-dotn/timestep);
			break;
		case MASS_MOMENTUM_CALL:
		case UPDATE_STRAINS_LAST_CALL:
			dotn = DotVectors(&pk, norm);
			AddPkScaled(norm, -dotn);
			break;
		case GRID_FORCES_CALL:
		{	// Start adjusting force for nodal velocity BC
			//      Ftot = (Ftot.norm) norm + (Ftot.tang) tang
			// But, now we want it to start with
			//      Ftotnew = -(pk.norm)/deltime norm + (Ftot.tang) tang
			//      Ftotnew = Ftot - ((Ftot.norm) + (pk.norm)/deltime) norm
			// Note that if >1 BC with same norm on same node, the net result after first pass is
			//      Ftot1 = -(pk.norm)/deltime norm + (Ftot.tang) tang
			// Then on second pass, change in force is:
			//		dotf = Ftot1.norm = -(pk.norm)/deltime
			//      -dotf-dotp/deltime = (pk.norm)/deltime - (pk.norm)/deltime = 0
			// But might have physical issues if components of norm overlap on the same node such as
			//    	x axis and  skew xy or xz on the same node
			double dotf = DotVectors(&ftot, norm);
			double dotp = DotVectors(&pk, norm);
			// the change in force is (-dotf-dotp/deltime)*norm, which is zero for second BC with same norm
			Vector deltaF;
			CopyScaleVector(&deltaF,norm,-dotf-dotp/deltime);
			AddFtot(&deltaF);
			AddVector(freaction,&deltaF);
			break;
		}
		case UPDATE_GRID_STRAINS_CALL:
			dotn = DotVectors(vk, norm);
			AddScaledVector(vk, norm, -dotn);
			break;
		default:
			break;
	}
}

// add one component of velocity (FMPM only), momentum, or force
void MatVelocityField::AddVelocityBC(Vector *norm,double vel,int passType,double deltime,Vector *freaction)
{
	switch(passType)
	{	case MASS_MOMENTUM_CALL:
			AddPkScaled(norm, mass*vel);
// #if ADJUST_COPIED_PK == 2
			// vk[pkCopy] = pk;
// #endif
			break;
		case UPDATE_MOMENTUM_CALL:
		{	double pvel = mass*vel;
			AddPkScaled(norm,pvel);
			AddFtotScaled(norm,pvel/timestep);
			break;
		}
		case UPDATE_STRAINS_LAST_CALL:
			AddPkScaled(norm, mass*vel);
			break;
		case GRID_FORCES_CALL:
		{	// the change in force is (mass*vel/deltime)*norm
			Vector deltaF;
			CopyScaleVector(&deltaF,norm,mass*vel/deltime);
			AddFtot(&deltaF);
			AddVector(freaction,&deltaF);
			break;
		}
		case UPDATE_GRID_STRAINS_CALL:
			AddScaledVector(vk, norm, vel);
			break;
		default:
			break;
	}
}

//#if ADJUST_COPIED_PK==1
// set symmetry plane momenta to zero in copied momenta
void MatVelocityField::AdjustForSymmetryBC(int fixedDirection)
{   if(fixedDirection&XSYMMETRYPLANE_DIRECTION)
 		vk[pkCopy].x=0.;
    if(fixedDirection&YSYMMETRYPLANE_DIRECTION)
		vk[pkCopy].y=0.;
    if(fixedDirection&ZSYMMETRYPLANE_DIRECTION)
		vk[pkCopy].z=0.;
}
//#endif

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
{   if(mvf==NULL) return false;										// no field
	if(mvf->IsRigidField()) return false;							// it is rigid
	if(fieldNum>0 && mvf->IgnoresCracks()) return false;			// field is not in memory
	return mvf->numberPoints>0;										// true if active
}

// return true if references field that is active, is not rigid, and matches requireCracks
// (i.e., if requireCracks is true, this field must see cracks, otherwise any field is OK)
// In NairnMPM, same as having any nonrigid particles
bool MatVelocityField::ActiveNonrigidSeesCracksField(MatVelocityField *mvf,bool requireCracks)
{   if(mvf==NULL) return false;										// no field
	if(mvf->numberPoints==0 || mvf->IsRigidField()) return false;	// no points or rigid
	return requireCracks ? !mvf->IgnoresCracks() : true ;
}
