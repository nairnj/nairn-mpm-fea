/********************************************************************************
	CrackVelocityField.cpp
	nairn-mpm-fea

	Created by John Nairn on 11 August 2009.
	Copyright (c) 2009 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Nodes/CrackVelocityField.hpp"
#include "Nodes/CrackVelocityFieldSingle.hpp"
#include "Nodes/CrackVelocityFieldMulti.hpp"
#include "Exceptions/CommonException.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Materials/MaterialBase.hpp"
#include "Materials/RigidMaterial.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Global_Quantities/BodyForce.hpp"

#pragma mark INITIALIZATION

// Constructors - vfld is field ID and cnum is crack number
// throws std::bad_alloc
CrackVelocityField::CrackVelocityField(int num,short theLoc,int cnum)
{
	// field number [0] to [3]
	fieldNum = num;
	
	// pointers for material velocity fields (1 if single material modes, or number of materials in use)
	mvf = new MatVelocityField *[maxMaterialFields];
	
	// set all to NULL, they are created as needed in task 1
	int i;
	for(i=0;i<maxMaterialFields;i++)
		mvf[i]=NULL;
	
	// in single material mode, always create the one material velocity field
	// (and it is never rigid and never ignores cracks for flags is zero
	if(maxMaterialFields==1)
	{	mvf[0] = CreateMatVelocityField(0);
	}

	// clear all data.
    // But creating it, implies it has points so set hasCrackPoints to true
	Zero(theLoc,cnum,false);
	hasCrackPoints = true;
	df = NULL;
	cm = NULL;
}

// Destructor - child class deletes material velocity fields as needed
CrackVelocityField::~CrackVelocityField()
{	
	delete [] mvf;
	DeleteStrainField();
}

// Create the type of material velocity field that is needed
// throws std::bad_alloc
MatVelocityField *CrackVelocityField::CreateMatVelocityField(int fieldFlags)
{
	return new MatVelocityField(fieldFlags);
}

#pragma mark TASK 0 METHODS

// zero this velocity field, set loc and crackNum for first crack and optionally
//		zero the material velocity fields
// WARNING: if zeroMVFS is true (which zeros velocity field too), must
//      NULL out fields that ignore crack before calling this method
void CrackVelocityField::Zero(short theLoc,int cnum,bool zeroMVFs)
{
	// duplicate this part in contructor
	loc[FIRST_CRACK]=theLoc;
	loc[SECOND_CRACK]=NO_CRACK;
	crackNum[FIRST_CRACK]=cnum;
	crackNum[SECOND_CRACK]=0;
	ZeroVector(&norm[FIRST_CRACK]);
	ZeroVector(&norm[SECOND_CRACK]);
	numberPoints=0;
    hasCrackPoints = false;

	// can't call when constructing because subclass not ready yet for this
	// pure virtual method
	if(zeroMVFs) ZeroMatFields();
}

// Called in intitation to preallocate material velocituy fields
// throws std::bad_alloc
void CrackVelocityField::AddMatVelocityField(int matfld) {}

bool CrackVelocityField::NeedsMatVelocityField(int matfld) const { return FALSE; }

// Make sure this crack velocity field on ghost node matches a real one
// Called during time step initialization if needed (i.e., cracks and/or multimaterial mode)
// throws std::bad::alloc
void CrackVelocityField::MatchRealFields(CrackVelocityField *rcvf)
{	// crack details
	loc[0] = rcvf->loc[0];
	loc[1] = rcvf->loc[1];
	crackNum[0] = rcvf->crackNum[0];
	crackNum[1] = rcvf->crackNum[1];
	CopyVector(&norm[0],&rcvf->norm[0]);
	CopyVector(&norm[1],&rcvf->norm[1]);
	MatchMatVelocityFields(rcvf->GetMaterialVelocityFields());
}

// match material velocity fields (only needed in multimaterial node
void CrackVelocityField::MatchMatVelocityFields(MatVelocityField **rmvf) {}

#pragma mark TASK 1 and 6 METHODS

// Called only in task 1 - add to momemtum and count number of material points
// If first material point for this velocity field, save the velocity too
void CrackVelocityField::AddMomentumTask1(int matfld,Vector *addPk,Vector *vel,int numPts)
{	// should have been allocated if needed in initialization step
	mvf[matfld]->AddMomentumTask1(addPk,vel,numPts);
	
	// count points in this crack velocity field during task 1
	numberPoints += numPts;
}

// Called only in task 1 - add rigid velcity
// Temporarily uses mvf[0]->pk, mvf[0]->vk[0], numberPoints
// This only used when extrapolating rigid BCs
void CrackVelocityField::AddRigidVelocityAndFlags(Vector *addPk,double fnmp,int setFlags)
{
	// save momentum if any direction is controlled
	if(setFlags & CONTROL_ANY_DIRECTION)
	{	// here addPk = fnmp*vel
		AddVector(&mvf[0]->pk,addPk);
	
		// add separate mass for each velocity that is set
		Vector *rmass = mvf[0]->vk;
		if(setFlags & CONTROL_X_DIRECTION) rmass->x += fnmp;
		if(setFlags & CONTROL_Y_DIRECTION) rmass->y += fnmp;
		if(setFlags & CONTROL_Z_DIRECTION) rmass->z += fnmp;
	}
	
	// collect flags in nummberPoints
	numberPoints |= setFlags;
}

// Read extrapolation get the velocity
// Restore mvf[0]->pk, mvf[0]->vk[0], and numberPoints to zero
int CrackVelocityField::ReadAndZeroRigidVelocity(Vector *rvel)
{
	// if not used, return 0
	if(numberPoints==0) return 0;
	
	// get flags and set any controlled velocities
	int tempFlags = numberPoints;
	if(tempFlags & CONTROL_ANY_DIRECTION)
	{	Vector *rpk = &mvf[0]->pk;
		Vector *rmass = mvf[0]->vk;
		if(tempFlags & CONTROL_X_DIRECTION) rvel->x = rpk->x/rmass->x;
		if(tempFlags & CONTROL_Y_DIRECTION) rvel->y = rpk->y/rmass->y;
		if(tempFlags & CONTROL_Z_DIRECTION) rvel->z = rpk->z/rmass->z;
		
		// zero used data
		ZeroVector(rpk);
		ZeroVector(rmass);
	}
	
	// return and return set directions
	numberPoints = 0;
	return tempFlags;
}

// add to mass (task 1) and field was allocated (if needed) in AddCrackVelocityField()
void CrackVelocityField::AddMass(int matfld,double mnode) { mvf[matfld]->mass += mnode; }

// add "mass" for  rigid particle (task 1) (only functions in CrackVelocityFieldMulti)
void CrackVelocityField::AddMassTask1(int matfld,double mnode,int numPts) { }

// Add to mass gradient (overridden in CrackVelocityFieldMulti where it is needed)
void CrackVelocityField::AddVolumeGradient(int matfld,MPMBase *mptr,double dNdx,double dNdy,double dNdz) {}
void CrackVelocityField::AddVolumeGradient(int matfld,Vector *grad) {}

// Copy mass and momentum from ghost to real node
void CrackVelocityField::CopyMassAndMomentum(NodalPoint *real)
{	if(mvf[0]!=NULL) mvf[0]->CopyMassAndMomentum(real,fieldNum,0);
}

// Add to momentum vector (second pass after fields so not need to count points)
void CrackVelocityField::AddMomentumTask6(int matfld,double wt,Vector *vel)
{	// momentum
	AddScaledVector(&mvf[matfld]->pk,vel,wt);		// in g mm/sec
}

// Copy mass and momentum from ghost to real node
void CrackVelocityField::CopyMassAndMomentumLast(NodalPoint *real)
{	if(mvf[0]!=NULL) mvf[0]->CopyMassAndMomentumLast(real,fieldNum,0);
}

#pragma mark TASK 3 METHODS

// Add to internal force
void CrackVelocityField::AddFtotTask3(int matfld,Vector *f)
{	mvf[matfld]->AddFtot(f);
}

// Copy grid forces from ghost node to real node (nonrigid only)
void CrackVelocityField::CopyGridForces(NodalPoint *real)
{	if(mvf[0]!=NULL) mvf[0]->CopyGridForces(real,fieldNum,0);
}

#pragma mark TASK 5 METHODS

// add to vStarNext
void CrackVelocityField::AddVStarNext(int matfld,Vector *vStarPrevj,double weight,double weightContact)
{	mvf[matfld]->AddVStarNext(vStarPrevj,weight,weightContact);
}

// Get vStarPrev pointer
Vector *CrackVelocityField::GetVStarPrev(int matfld) const
{	return mvf[matfld]->GetVStarPrev();
}

// Get real node mass for material
double CrackVelocityField::GetMaterialMass(int matfld) const
{	return mvf[matfld]->mass;
}

// Increment velocity and acceleration for this material point using one velocity field which must be there
void CrackVelocityField::IncrementDelvaTask5(int matfld,double fi,GridToParticleExtrap *gp) const
{   mvf[matfld]->IncrementNodalVelAcc(fi,gp);
}

#pragma mark TASK 6 METHODS

#pragma mark TASK 7 J AND K CALCULATION METHODS

// create strain field for this crack velocity field
// throws std::bad_alloc
void CrackVelocityField::CreateStrainField(void)
{	df=new DispField;
	
	df->du.x=0.;
	df->du.y=0.;
	df->dv.x=0.;
	df->dv.y=0.;
	df->kinetic=0.;
	df->work=0.;
	df->stress.xx=0.;
	df->stress.yy=0.;
	df->stress.zz=0.;
	df->stress.xy=0.;
    df->mass=0.;
    df->vx=0.;
    df->vy=0.;
	
	// if more than one material, track material type to set crack tip material changes
	if(numActiveMaterials>1)
	{	df->matWeight = new double[numActiveMaterials];
		int i;
		for(i=0;i<numActiveMaterials;i++)
			df->matWeight[i] = 0.;
	}
	else
		df->matWeight=NULL;
}

// create strain field for this crack velocity field
void CrackVelocityField::DeleteStrainField(void)
{	if(df!=NULL)
	{	if(df->matWeight!=NULL) delete [] df->matWeight;
		delete df;
		df=NULL;
	}
}

#pragma mark TASK 8 METHODS

// Increment velocity when moving crack surfaces.
// This called when moving one side of the crack in velocity field for that
//   side of the crack. Currently looks only at materials that see the cracks
//   which is analgous to way it works for rigid material inside a crack
// On output, *fieldMass is fi*total mass to allow mass-weighted extrapolations
//   Older code did not use mass weighted but that appeared less stable if a node
//   had very low mass and inaccurate velocity
// Older code optionally screened out low mass nodes, but that was dropped in favor of
//   always doing mass weighted sums
short CrackVelocityField::IncrementDelvTask8(double fi,Vector *delV,Vector *dela,double *fieldMass)
{
	// get CM momentum, force, and total mass this crack velocity field
    // In multimaterial mode, it combines all materials on this side of the crack
    // For single crack, just gets momentum, force, and mass of that field
	double totalMass;
	bool hasParticles;
	Vector totalFtot;
	ZeroVector(&totalFtot);
	Vector totalPk = GetCMatMomentum(hasParticles,&totalMass,&totalFtot,false);

	// skip no particles
	if(!hasParticles) return false;
	
	// increment velocity extrapolation by Sip pi
	AddScaledVector(delV,&totalPk,fi);

	// increment acceleration by Sip fi
	AddScaledVector(dela,&totalFtot,fi);
	
	// mass weighted normalization
	*fieldMass = fi*totalMass;

	return true;
}

// Collect momenta and add to vector when finding CM velocity to move crack planes
// Also increment the mass
// return if found nonrigid points that see cracks
bool CrackVelocityField::CollectMomentaInCrackField(Vector *totalPk,double *velocityMass,Vector *totalFtot,bool useVelocity) const
{	bool hasParticles;
	double fieldMass;
	Vector fieldFtot;
	ZeroVector(&fieldFtot);
	Vector fieldPk = GetCMatMomentum(hasParticles,&fieldMass,totalFtot,useVelocity);
	if(hasParticles)
	{	AddVector(totalPk,&fieldPk);
		if(totalFtot!=NULL) AddVector(totalFtot,&fieldFtot);
		*velocityMass += fieldMass;
	}
	return hasParticles;
}

// Collect momenta and add to vector when finding CM velocity to move crack planes
void CrackVelocityField::SetCMVelocityTask8(Vector *velCM,int totalParticles,Vector *accCM)
{	// store in normal vector and crack number to save memory (normal not needed again in this time step)
	if(cm==NULL) cm = new CenterMassField;
	cm->numParticles = totalParticles;
	cm->cmVel = *velCM;
	cm->cmAcc = *accCM;
}

// Return CM velocity for crack updates
bool CrackVelocityField::GetCMVelocityTask8(Vector *velCM,Vector *accCM) const
{	// stored in normal vector and crack number to save memory (they are not needed again in this time step)
	if(cm==NULL) return false;
	if(cm->numParticles==0) return false;
	*velCM = cm->cmVel;
	*accCM = cm->cmAcc;
	return true;
}

#pragma mark MATERIAL CONTACT AND INTERFACES IN SUBCLASSES

// Called in multimaterial mode to check contact at nodes with multiple materials
void CrackVelocityField::MaterialContactOnCVF(MaterialContactNode *mcn,double deltime,int callType)
{	
}

// retrieve mass gradient (overridden in CrackVelocityFieldMulti where it is needed
void CrackVelocityField::GetVolumeGradient(int matfld,const NodalPoint *ndptr,Vector *grad,double scale) const { ZeroVector(grad); }

// Adjust vector for symmetry planes, if keepNormalized, renormalize on any change
// Component of vector normal to symmetry plane direction is set to zero
void CrackVelocityField::AdjustForSymmetry(NodalPoint *ndptr,Vector *norm,bool keepNormalized) const
{
    // if has any, have to check each one
    bool renormalize = false;
    if(ndptr->fixedDirection&XSYMMETRYPLANE_DIRECTION)
    {   norm->x = 0.;
        renormalize = true;
    }
    if(ndptr->fixedDirection&YSYMMETRYPLANE_DIRECTION)
    {   norm->y = 0.;
        renormalize = true;
    }
    if(ndptr->fixedDirection&ZSYMMETRYPLANE_DIRECTION)
    {   norm->z = 0.;
        renormalize = true;
    }
    if(renormalize && keepNormalized)
        CopyScaleVector(norm,norm,1./sqrt(DotVectors(norm,norm)));
}

#pragma mark PROPERTIES FOR CRACK AND MATERIAL CONTACT

// add to normal vector
void CrackVelocityField::AddNormals(Vector *cnorm,int which)
{	
	AddVector(&norm[which], cnorm);
}

// Add displacements and position (getting different as needed)
void CrackVelocityField::AddVolumeDisplacement(int matfld,double wtVol,double wt,Vector ppos,Vector *pdisp)
{	mvf[matfld]->AddContactVolume(wtVol);
	if(mpmgrid.positionIndex>=0)
		mvf[matfld]->AddContactVector(mpmgrid.positionIndex,&ppos,wt);
	if(mpmgrid.displacementIndex>=0)
		mvf[matfld]->AddContactVector(mpmgrid.displacementIndex,SubVector(&ppos,pdisp),wt);
}

// Add contact terms to selected field
void CrackVelocityField::AddContactTerms(int matfld,ContactTerms *contactInfo)
{	if(mpmgrid.positionIndex>=0)
		mvf[matfld]->AddContactVector(mpmgrid.positionIndex,&contactInfo->terms[mpmgrid.positionIndex]);
	if(mpmgrid.displacementIndex>=0)
		mvf[matfld]->AddContactVector(mpmgrid.displacementIndex,&contactInfo->terms[mpmgrid.displacementIndex]);
	mvf[matfld]->AddContactVolume(contactInfo->cvolume);
	if(mpmgrid.volumeGradientIndex>=0)
		mvf[matfld]->AddContactVector(mpmgrid.volumeGradientIndex,&contactInfo->terms[mpmgrid.volumeGradientIndex]);
}

#pragma mark ACCESSORS

// location for crack in this field (do not call for field [0] and may be undefined)
short CrackVelocityField::location(int which) { return loc[which]; }

// crack number for this field (do not call for field [0] and may be undefined)
int CrackVelocityField::crackNumber(int which) { return crackNum[which]; }

// if has crack matching supplied number and side, return the other crack, else return -1
// only make sense to call for field [3]
int CrackVelocityField::OppositeCrackTo(int cnum,int side,int *otherSide)
{	if(crackNum[FIRST_CRACK]==cnum)
	{	if(side==loc[FIRST_CRACK])
		{	*otherSide = loc[SECOND_CRACK];
			return crackNum[SECOND_CRACK];
		}
	}
	else if(crackNum[SECOND_CRACK]==cnum)
	{	if(side==loc[SECOND_CRACK])
		{	*otherSide = loc[FIRST_CRACK];
			return crackNum[FIRST_CRACK];
		}
	}
	return -1;
}
		
// set both at once
void CrackVelocityField::SetLocationAndCrack(short vfld,int cnum,int which)
{	loc[which]=vfld;
	crackNum[which]=cnum;
    hasCrackPoints = true;
}

// Get material velocity fields
MatVelocityField **CrackVelocityField::GetMaterialVelocityFields(void) { return mvf; }

// Get material velocity fields
MatVelocityField *CrackVelocityField::GetMaterialVelocityField(int matfld) { return mvf[matfld]; }

// Get velocity for selected material field
Vector CrackVelocityField::GetVelocity(int matfld)
{	return mvf[matfld]->GetVelocity();
}

// total number of points (rigid and nonrigid included and mirrored are included)
// These apply to single material mode, which cannot have rigid materials
int CrackVelocityField::GetNumberPoints(void) { return numberPoints; }
void CrackVelocityField::SetNumberPoints(int numpts) { numberPoints = numpts; }
int CrackVelocityField::HasPointsThatSeeCracks(void) { return numberPoints>0 ? 1 : 0 ; }
int CrackVelocityField::GetNumberNonrigidMaterials(void) { return numberPoints>0 ? 1 : 0 ; }
int CrackVelocityField::GetNumberMaterials(void) { return 1; }

// Look for presence of nonrigid poitns (override in CrackVelocityFieldMulti)
bool CrackVelocityField::HasPointsNonrigid(void) const { return numberPoints>0; }

// for debugging
void CrackVelocityField::Describe(void) const
{
	cout << "# Crack Field: npts="<<  numberPoints << " mass=" << GetTotalMass(true) << " cracking mass=" << GetTotalMass(false)
		<< " vol=(" << GetContactVolumeNonrigid(false) << "," << GetContactVolumeNonrigid(false) << ")" << endl;
	if(crackNum[0]>0)
	{	cout << "#     crack 1=#" << crackNum[0] << ", loc=";
		if(loc[0]==ABOVE_CRACK) cout << "above"; else cout << "below";
		PrintVector(", n=",&norm[0]);
		cout << endl;
	}
	else
		cout << "#     non-crossing crack field" << endl;
	if(crackNum[1]>0)
	{	cout << "#     crack 2=#" << crackNum[1] << ", loc=" << loc[1];
		if(loc[0]==ABOVE_CRACK) cout << "above"; else cout << "below";
		PrintVector(", n=",&norm[1]);
		cout << endl;
	}
}

// add contact force on rigid material to the input vector
void CrackVelocityField::SumAndClearRigidContactForces(Vector *fcontact,bool clearForces,double scale,Vector *ftotal) {}

// return field number [0] to [3]
int CrackVelocityField::GetFieldNum(void) const { return fieldNum; }

#pragma mark CLASS METHODS

// return true if referenced field is active in this time step
bool CrackVelocityField::ActiveField(CrackVelocityField *cvf)
{ return cvf==NULL ? false : (cvf->numberPoints>0) ; }

// return true if this crack velocity field has material field active
bool CrackVelocityField::HasActiveMatField(CrackVelocityField *cvf,int matfld)
{	if(cvf==NULL) return false;
	if(cvf->numberPoints<=0) return false;
	MatVelocityField *mvf = cvf->GetMaterialVelocityField(matfld);
	return MatVelocityField::ActiveNonrigidField(mvf);
}

// return true if referenced field is active in this time step during velocity field allocation
bool CrackVelocityField::ActiveCrackField(CrackVelocityField *cvf)
{ return cvf==NULL ? false : cvf->hasCrackPoints ; }

// return true if referenced field is active AND has some non-rigid particles
bool CrackVelocityField::ActiveNonrigidField(CrackVelocityField *cvf)
{ return cvf==NULL ? false : cvf->HasPointsNonrigid() ; }

// return true if referenced field is active AND has some non-rigid particles AND first crack matches crack number
// do not call with cvf[0] because that will not track crack number
bool CrackVelocityField::ActiveNonrigidField(CrackVelocityField *cvf,int number)
{	if(cvf==NULL) return false;
	if(cvf->HasPointsNonrigid())
	{	if(cvf->crackNum[FIRST_CRACK]==number)
			return true;
	}
	return false;
}

// Single material mode gets CrackVelocityFieldSingle
// Multimaterial mode gets CrackVelocityFieldMulti
// CrackVelocityField is an abstract class
// throws std::bad_alloc
CrackVelocityField *CrackVelocityField::CreateCrackVelocityField(int fieldNum,short theLoc,int cnum)
{	if(maxMaterialFields==1)
		return (CrackVelocityField *)(new CrackVelocityFieldSingle(fieldNum,theLoc,cnum));
	else
		return (CrackVelocityField *)(new CrackVelocityFieldMulti(fieldNum,theLoc,cnum));
}

