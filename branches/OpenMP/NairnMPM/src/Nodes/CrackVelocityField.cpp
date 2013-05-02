/********************************************************************************
	CrackVelocityField.cpp
	NairnMPM

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
CrackVelocityField::CrackVelocityField(short theLoc,int cnum)
{
	// pointers for material velocity fields (1 if single material modes, or number of materials in use)
	mvf=(MatVelocityField **)malloc(sizeof(MatVelocityField *)*maxMaterialFields);
	if(mvf==NULL) throw CommonException("Memory error allocating material velocity field pointers.",
                                        "CrackVelocityField::CrackVelocityField");
	// set all to NULL, they are created as needed in task 1
	int i;
	for(i=0;i<maxMaterialFields;i++)
		mvf[i]=NULL;
	
	// in single material mode, always create the one material velocity field
	if(maxMaterialFields==1)
	{	mvf[0]=new MatVelocityField(FALSE);
		if(mvf[0]==NULL) throw CommonException("Memory error allocating material velocity field.",
											   "CrackVelocityField::CrackVelocityField");
	}
		
	// clear all data (except hasCrackPoints)
	Zero(theLoc,cnum,FALSE);
	hasCrackPoints=TRUE;
	df=NULL;
}

// Destructor
CrackVelocityField::~CrackVelocityField()
{	// delete any allocated material velocity fields
	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(mvf[i]!=NULL)
			delete mvf[i];
	}
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
{	if(mvf[matfld]==NULL) throw "encoutered missing material velocity field";
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
void CrackVelocityField::CopyMassAndMomentum(NodalPoint *real,int vfld)
{	if(mvf[0]!=NULL) mvf[0]->CopyMassAndMomentum(real,vfld,0);
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
void CrackVelocityField::CopyMassAndMomentumLast(NodalPoint *real,int vfld)
{	if(mvf[0]!=NULL) mvf[0]->CopyMassAndMomentumLast(real,vfld,0);
}

#pragma mark TASK 3 METHODS

// Add to internal force
void CrackVelocityField::AddFtotTask3(int matfld,Vector *f) { mvf[matfld]->AddFtot(f); }

// Copy grid forces from ghost node to real node (nonrigid only)
void CrackVelocityField::CopyGridForces(NodalPoint *real,int vfld)
{	if(mvf[0]!=NULL) mvf[0]->CopyGridForces(real,vfld,0);
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
	
	// if more than one material, track material tip to set crack tip material chnages
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

// Increment velocity when moving crack surface. This called after CM and total mass
// is stored in the first non-empty material field
short CrackVelocityField::IncrementDelvTask8(double fi,Vector *delV,double *fieldMass)
{	// skip low mass node
	double totalMass=GetTotalMass();
	if(totalMass/(*fieldMass)<1.e-6) return FALSE;			// skip low mass
	Vector totalPk=GetCMatMomentum();
	AddScaledVector(delV,&totalPk,fi/totalMass);			// increment
	*fieldMass=totalMass;
	return TRUE;
	/*
	*fieldMass=GetTotalMass();
	Vector totalPk=GetCMatMomentum();
	AddScaledVector(delV,&totalPk,fi);			// increment
	return TRUE;
	*/
}

// Collect momenta and add to vector when finding CM velocity to move crack planes
// return number of nonrigid points
int CrackVelocityField::CollectMomentaTask8(Vector *totalPk)
{	Vector fieldPk=GetCMatMomentum();
	AddVector(totalPk,&fieldPk);
	return GetNumberPointsNonrigid();
}

// Collect momenta and add to vector when finding CM velocity to move crack planes
void CrackVelocityField::SetCMVelocityTask8(Vector *velCM,int totalParticles)
{	// store in normal vector and crack number to save memory (normal not needed again in this time step)
	norm[0]=*velCM;
	crackNum[0]=totalParticles;
}

// Return CM velocity for crack updates
bool CrackVelocityField::GetCMVelocityTask8(Vector *velCM)
{	// stored in normal vector and crack number to save memory (they are not needed again in this time step)
	if(crackNum[0]==0) return FALSE;
	*velCM=norm[0];
	return TRUE;
}

#pragma mark MATERIAL CONTACT AND INTERFACES IN SUBCLASSES

// Called in multimaterial mode to check contact at nodes with multiple materials
void CrackVelocityField::MaterialContactOnCVF(NodalPoint *ndptr,int vfld,double deltime,bool postUpdate,
											  MaterialInterfaceNode **first,MaterialInterfaceNode **last)
{ return;
}

// retrieve mass gradient (overridden in CrackVelocityFieldMulti where it is needed
void CrackVelocityField::GetVolumeGradient(int matfld,const NodalPoint *ndptr,Vector *grad,double scale) const { ZeroVector(grad); }

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

// if has crack matching supplied number and side, return the other crack, else return -1
int CrackVelocityField::OppositeCrackTo(int cnum,int side)
{	if(crackNum[FIRST_CRACK]==cnum)
		return (side==loc[FIRST_CRACK]) ? crackNum[SECOND_CRACK] : -1 ;
	else if(crackNum[SECOND_CRACK]==cnum)
		return (side==loc[SECOND_CRACK]) ? crackNum[FIRST_CRACK] : -1 ;
	else
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

// Get velocity for selected material field
Vector CrackVelocityField::GetVelocity(int matfld)
{	return mvf[matfld]->GetVelocity();
}

// Get velocity for selected material field
Vector CrackVelocityField::GetContactForce(int matfld)
{	return mvf[matfld]->GetFtot();
}

// total number of points (rigid and nonrigid included)
int CrackVelocityField::GetNumberPoints(void) { return numberPoints; }

// total number of non-rigid points (override in CrackVelocityFieldMulti)
int CrackVelocityField::GetNumberPointsNonrigid(void) { return numberPoints; }

// for debugging
void CrackVelocityField::Describe(void)
{
	cout << "# Crack Field: npts="<<  numberPoints << " mass=" << GetTotalMass()
		<< " vol=" << GetVolumeTotal(1.) << endl;
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
void CrackVelocityField::SumAndClearRigidContactForces(Vector *fcontact,bool) {}

#pragma mark CLASS METHODS

// return true if referenced field is active in this time step
bool CrackVelocityField::ActiveField(CrackVelocityField *cvf)
{ return cvf==NULL ? (bool)FALSE : (cvf->numberPoints>0) ; }

// return true if referenced field is active in this time step during velocity field allocation
bool CrackVelocityField::ActiveCrackField(CrackVelocityField *cvf)
{ return cvf==NULL ? (bool)FALSE : cvf->hasCrackPoints ; }

// return true if referenced field is active AND has some non-rigid particles
bool CrackVelocityField::ActiveNonrigidField(CrackVelocityField *cvf)
{ return cvf==NULL ? (bool)FALSE : (cvf->GetNumberPointsNonrigid()>0) ; }

// create single or multi material crack velocity field as needed
CrackVelocityField *CrackVelocityField::CreateCrackVelocityField(short theLoc,int cnum)
{	if(maxMaterialFields==1)
		return (CrackVelocityField *)(new CrackVelocityFieldSingle(theLoc,cnum));
	else
		return (CrackVelocityField *)(new CrackVelocityFieldMulti(theLoc,cnum));
}




