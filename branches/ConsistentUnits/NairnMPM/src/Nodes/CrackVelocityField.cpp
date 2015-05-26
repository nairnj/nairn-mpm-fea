/********************************************************************************
	CrackVelocityField.cpp
	nairn-mpm-fea

	Created by John Nairn on 11 August 2009.
	Copyright (c) 2009 John A. Nairn, All rights reserved.
********************************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Nodes/CrackVelocityField.hpp"
#include "Nodes/CrackVelocityFieldSingle.hpp"
#include "Nodes/CrackVelocityFieldMulti.hpp"
#include "Exceptions/CommonException.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Materials/MaterialBase.hpp"

#pragma mark INITIALIZATION

// Constructors - vfld is field ID and cnum is crack number
CrackVelocityField::CrackVelocityField(int num,short theLoc,int cnum)
{
	// field number [0] to [3]
	fieldNum = num;
	
	// pointers for material velocity fields (1 if single material modes, or number of materials in use)
	mvf=(MatVelocityField **)malloc(sizeof(MatVelocityField *)*maxMaterialFields);
	if(mvf==NULL) throw CommonException("Memory error allocating material velocity field pointers.",
                                        "CrackVelocityField::CrackVelocityField");
	// set all to NULL, they are created as needed in task 1
	int i;
	for(i=0;i<maxMaterialFields;i++)
		mvf[i]=NULL;
	
	// in single material mode, always create the one material velocity field
	// (and it is never rigid and never ignores cracks for flags is zero
	if(maxMaterialFields==1)
	{	mvf[0]=new MatVelocityField(0);
		if(mvf[0]==NULL) throw CommonException("Memory error allocating material velocity field.",
											   "CrackVelocityField::CrackVelocityField");
	}
		
	// clear all data.
    // But creating it, implies it has points so set hasCrackPoints to TRUE
	Zero(theLoc,cnum,FALSE);
	hasCrackPoints=TRUE;
	df=NULL;
}

// Destructor - child class deletes material velocity fields as needed
CrackVelocityField::~CrackVelocityField()
{	
	free(mvf);
	DeleteStrainField();
}

#pragma mark TASK 0 METHODS

// zero this velocity field, set loc and crackNum for first crack and optionally
//		zero or delete the material velocity fields
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
    hasCrackPoints=FALSE;
	
	// can't call when constructing because subclass not ready yet for this
	// pure virtual method
	if(zeroMVFs) ZeroMatFields();
}

// Called in intitation to preallocate material velocituy fields
void CrackVelocityField::AddMatVelocityField(int matfld) {}
bool CrackVelocityField::NeedsMatVelocityField(int matfld) const { return FALSE; }

// Make sure this crack velocity field on ghost node matches a real one
// Called during time step initialization if needed (i.e., cracks and/or multimaterial mode)
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

// add to mass (task 1) and field was allocated (if needed) in AddCrackVelocityField()
void CrackVelocityField::AddMass(int matfld,double mnode) { mvf[matfld]->mass += mnode; }

// add "mass" for  rigid particle (task 1) (only functions in CrackVelocityFieldMulti)
void CrackVelocityField::AddMassTask1(int matfld,double mnode,int numPts) { }

// Add to mass gradient (overridden in CrackVelocityFieldMulti where it is needed)
void CrackVelocityField::AddVolumeGradient(int matfld,MPMBase *mptr,double dNdx,double dNdy,double dNdz) {}

// Add to mass gradient (overridden in CrackVelocityFieldMulti where it is needed)
void CrackVelocityField::CopyVolumeGradient(int matfld,Vector *grad) {}

// Copy mass and momentum from ghost to real node
void CrackVelocityField::CopyMassAndMomentum(NodalPoint *real)
{	if(mvf[0]!=NULL) mvf[0]->CopyMassAndMomentum(real,fieldNum,0);
}

// Add to momentum vector (second pass after fields so not need to count points)
void CrackVelocityField::AddMomentumTask6(int matfld,double wt,Vector *vel)
{	// momentum
	AddScaledVector(&mvf[matfld]->pk,vel,wt);		// in g mm/sec
    
    // save velocity if only one point on this node
    if(numberPoints==1)
        mvf[matfld]->SetVelocity(vel);
}

// Copy mass and momentum from ghost to real node
void CrackVelocityField::CopyMassAndMomentumLast(NodalPoint *real)
{	if(mvf[0]!=NULL) mvf[0]->CopyMassAndMomentumLast(real,fieldNum,0);
}

#pragma mark TASK 3 METHODS

// Add to internal force
void CrackVelocityField::AddFtotTask3(int matfld,Vector *f) { mvf[matfld]->AddFtot(f); }

// Copy grid forces from ghost node to real node (nonrigid only)
void CrackVelocityField::CopyGridForces(NodalPoint *real)
{	if(mvf[0]!=NULL) mvf[0]->CopyGridForces(real,fieldNum,0);
}

#pragma mark TASK 5 METHODS

// Increment velocity and acceleration for this material point using one velocity field which must be there
void CrackVelocityField::IncrementDelvaTask5(int matfld,double fi,Vector *delv,Vector *dela) const
{   mvf[matfld]->IncrementNodalVelAcc(fi,delv,dela);
}

#pragma mark TASK 6 METHODS

#pragma mark TASK 7 J AND K CALCULATION METHODS

// create strain field for this crack velocity field
void CrackVelocityField::CreateStrainField(void)
{	df=new DispField;
    if(df==NULL) throw CommonException("Memory error allocating new strain field.",
										"CrackVelocityField::CreateStrainField");
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
	
	// if more than one material, track material tip to set crack tip material changes
	if(numActiveMaterials>1)
	{	df->matWeight=(double *)malloc(sizeof(double)*numActiveMaterials);
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
	{	if(df->matWeight!=NULL) delete df->matWeight;
		delete df;
		df=NULL;
	}
}

#pragma mark TASK 8 METHODS

// Increment velocity when moving crack surface.
// This called when moving one side of the crack in velocity field for that
//   side of the crack. Currently looks only at materials that sees the cracks
//   which is analgous to way it works for rigid material inside a crack
short CrackVelocityField::IncrementDelvTask8(double fi,Vector *delV,double *fieldMass)
{	
#ifndef CRACK_SURFACE_BY_MOMENTUM_EXTRAP
	// skip if no mass on this node
	if(*fieldMass==0.) return false;
#endif
	
	// get CM momentum
	double totalMass;
	bool hasParticles;
	Vector totalPk = GetCMatMomentum(hasParticles,&totalMass);

#ifdef CRACK_SURFACE_BY_MOMENTUM_EXTRAP
	// skip no particles
	if(!hasParticles) return false;
    
    // increment momentum extrapolation by Sip pi
	AddScaledVector(delV,&totalPk,fi);
#else
	// skip no particles or low mass
	if(!hasParticles || totalMass/(*fieldMass)<1.e-5) return false;
    
    // increment velocity extrapolation by Sip vi
	AddScaledVector(delV,&totalPk,fi/totalMass);
#endif
	
	// return mass
	*fieldMass = totalMass;
	return true;
}

// Collect momenta and add to vector when finding CM velocity to move crack planes
// Also increment the mass
// return if found nonrigid points that see cracks
bool CrackVelocityField::CollectMomentaTask8(Vector *totalPk,double *velocityMass) const
{	bool hasParticles;
	double fieldMass;
	Vector fieldPk = GetCMatMomentum(hasParticles,&fieldMass);
	if(hasParticles)
	{	AddVector(totalPk,&fieldPk);
		*velocityMass += fieldMass;
	}
	return hasParticles;
}

// Collect momenta and add to vector when finding CM velocity to move crack planes
void CrackVelocityField::SetCMVelocityTask8(Vector *velCM,int totalParticles)
{	// store in normal vector and crack number to save memory (normal not needed again in this time step)
	norm[0]=*velCM;
	crackNum[0]=totalParticles;
}

// Return CM velocity for crack updates
bool CrackVelocityField::GetCMVelocityTask8(Vector *velCM) const
{	// stored in normal vector and crack number to save memory (they are not needed again in this time step)
	if(crackNum[0]==0) return false;
	*velCM=norm[0];
	return true;
}

#pragma mark MATERIAL CONTACT AND INTERFACES IN SUBCLASSES

// Called in multimaterial mode to check contact at nodes with multiple materials
void CrackVelocityField::MaterialContactOnCVF(NodalPoint *ndptr,double deltime,int callType,
											  MaterialInterfaceNode **first,MaterialInterfaceNode **last)
{ return;
}

// retrieve volume gradient for matnum (1 based) in crack field only (or zero if
// not there or not tracked (subclass overrides)
bool CrackVelocityField::HasVolumeGradient(int matfld) const { return FALSE; }

// retrieve mass gradient (overridden in CrackVelocityFieldMulti where it is needed
void CrackVelocityField::GetVolumeGradient(int matfld,const NodalPoint *ndptr,Vector *grad,double scale) const { ZeroVector(grad); }

// Adjust vector for symmetry planes, if keepNormalized, renormalize on any change
void CrackVelocityField::AdjustForSymmetry(NodalPoint *ndptr,Vector *norm,bool keepNormalized) const
{
    // if has any, have to check eachone
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

// add to normal vector (2D only for crack contact)
void CrackVelocityField::AddNormals(Vector *cnorm,int which)
{	norm[which].x+=cnorm->x;
	norm[which].y+=cnorm->y;
}

// Add displacements
void CrackVelocityField::AddDisplacement(int matfld,double wt,Vector *pdisp)
{	AddScaledVector(&mvf[matfld]->disp,pdisp,wt);
}

// Add volume
void CrackVelocityField::AddVolume(int matfld,double wtVol)
{	mvf[matfld]->AddContactVolume(wtVol);
}

#pragma mark ACCESSORS

// location for crack in this field
short CrackVelocityField::location(int which) { return loc[which]; }

// crack number for this field
int CrackVelocityField::crackNumber(int which) { return crackNum[which]; }

// location and crack for this field
bool CrackVelocityField::crackAndLocation(int which,int cnum,int side)
{	return crackNum[which]==cnum && loc[which]==side;
}

// if has crack matching supplied number and side, return the other crack, else return -1
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
    hasCrackPoints=TRUE;
}

// Get material velocity fields
MatVelocityField **CrackVelocityField::GetMaterialVelocityFields(void) { return mvf; }

// Get material velocity fields
MatVelocityField *CrackVelocityField::GetMaterialVelocityField(int matfld) { return mvf[matfld]; }

// Get velocity for selected material field
Vector CrackVelocityField::GetVelocity(int matfld)
{	return mvf[matfld]->GetVelocity();
}

// total number of points (rigid and nonrigid included)
int CrackVelocityField::GetNumberPoints(void) { return numberPoints; }
void CrackVelocityField::SetNumberPoints(int numpts) { numberPoints = numpts; }

// Look for presence of nonrigid poitns (override in CrackVelocityFieldMulti)
bool CrackVelocityField::HasPointsNonrigid(void) const { return numberPoints>0; }

// for debugging
void CrackVelocityField::Describe(void) const
{
	cout << "# Crack Field: npts="<<  numberPoints << " mass=" << GetTotalMass(true) << " cracking mass=" << GetTotalMass(false)
		<< " vol=" << GetVolumeTotal(NULL) << endl;
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

// return true if referenced field is active in this time step during velocity field allocation
bool CrackVelocityField::ActiveCrackField(CrackVelocityField *cvf)
{ return cvf==NULL ? false : cvf->hasCrackPoints ; }

// return true if referenced field is active AND has some non-rigid particles
bool CrackVelocityField::ActiveNonrigidField(CrackVelocityField *cvf)
{ return cvf==NULL ? false : cvf->HasPointsNonrigid() ; }

// return true if referenced field is active AND has some non-rigid particles AND first crack matches crack number
bool CrackVelocityField::ActiveNonrigidField(CrackVelocityField *cvf,int number)
{	if(cvf==NULL) return false;
	if(cvf->HasPointsNonrigid())
	{	if(cvf->crackNum[FIRST_CRACK]==number)
			return true;
	}
	return false;
}

// create single or multi material crack velocity field as needed
CrackVelocityField *CrackVelocityField::CreateCrackVelocityField(int fieldNum,short theLoc,int cnum)
{	if(maxMaterialFields==1)
		return (CrackVelocityField *)(new CrackVelocityFieldSingle(fieldNum,theLoc,cnum));
	else
		return (CrackVelocityField *)(new CrackVelocityFieldMulti(fieldNum,theLoc,cnum));
}




