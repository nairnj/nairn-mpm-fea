/********************************************************************************
    More NodalPoint.cpp for MPM only
    NairnMPM
    
    Created by John Nairn on Mar 17 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
	
	NodalPoint: Methods for MPM only
********************************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Nodes/NodalPoint.hpp"
#include "System/ArchiveData.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Exceptions/CommonException.hpp"
#include "Exceptions/MPMWarnings.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Cracks/CrackNode.hpp"
#include "Cracks/CrackSegment.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Nodes/CrackVelocityFieldMulti.hpp"
#include "Nodes/MatVelocityField.hpp"
#include "Nodes/MaterialInterfaceNode.hpp"

// class statics
double NodalPoint::interfaceEnergy=0.;

#pragma mark INITIALIZATION

// Destructor
NodalPoint::~NodalPoint()
{
	if(cvf!=NULL)
	{	int i;
		for(i=1;i<maxCrackFields;i++)
		{	if(cvf[i]!=NULL) delete cvf[i];
		}
		free(cvf);
	}
}

// Called by PreliminaryCalcs() just before the MPM analysis starts
// Can allocate things that were not known while reading the XML file
void NodalPoint::PrepareForFields()
{
	// need maxCrackFields Crack Velocity fields (1 if no cracks or MAX_FIELDS_FOR_CRACKS if any cracks)
	cvf=(CrackVelocityField **)malloc(sizeof(CrackVelocityField *)*maxCrackFields);
	if(cvf==NULL) throw CommonException("Memory error allocating crack velocity field pointers.",
										"NodalPoint::PrepareForFields");

	// always create first one, set rest (may be none) to NULL
	cvf[0]=CrackVelocityField::CreateCrackVelocityField(0,0);
	if(cvf[0]==NULL) throw CommonException("Memory error allocating crack velocity field 0.",
										   "NodalPoint::PrepareForFields");
	int i;
	for(i=1;i<maxCrackFields;i++) cvf[i]=NULL;
}


// zero data and reduce to one field at start of a step
void NodalPoint::InitializeForTimeStep(void)
{	
	for(int i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->Zero(0,0,TRUE);
	}
	gVolume=0.;
	gConcentration=0.;
	gTemperature=0.;
	gRhoVCp=0.;
	fdiff=0.;
	fcond=0.;
}

#pragma mark TASK 1 METHODS

// Add mass for selected field
void NodalPoint::AddMass(short vfld,int matfld,double mnode) { cvf[vfld]->AddMass(matfld,mnode); }

// for rigid particles, adding mass is counting number of rigid particles
void NodalPoint::AddMassTask1(short vfld,int matfld,double mnode) { cvf[vfld]->AddMassTask1(matfld,mnode); }

// Add to momentum vector (first pass - allocate cvf[] if needed) (both 2D and 3D)
short NodalPoint::AddMomentumTask1(int matfld,CrackField *cfld,double wt,Vector *vel)
{	short vfld=0;
	
	// only 1 or no cracks, with relevant settings in cfld[0]
	if(cfld[1].loc==NO_CRACK)
	{	switch(cfld[0].loc)
		{	case NO_CRACK:
				vfld=0;
				break;
			case ABOVE_CRACK:
			case BELOW_CRACK:
				// First try to use cvf[1]. It can be used if empty or if already there for same crack and location
				if(!CrackVelocityField::ActiveField(cvf[1]))
				{	// cvf[1] is empty
					if(cvf[1]==NULL)
					{	cvf[1]=CrackVelocityField::CreateCrackVelocityField(cfld[0].loc,cfld[0].crackNum);
						if(cvf[1]==NULL) throw CommonException("Memory error allocating crack velocity field 1.",
															   "NodalPoint::AddMomentumTask1");
					}
					else
						cvf[1]->SetLocationAndCrack(cfld[0].loc,cfld[0].crackNum,FIRST_CRACK);
					vfld=1;
				}
				
				// if the one crack crossed is same crack as cvf[1], then it should go in field [0] or field [1]
				else if(cvf[1]->crackNumber(FIRST_CRACK)==cfld[0].crackNum)
				{	if(cvf[1]->location(FIRST_CRACK)==cfld[0].loc)
					{	// found another point for [1]
						vfld=1;
					}
					else
					{	// it can only be field 0
						vfld=0;
						
						// Here means both above and below crack for field [1], which can only happen is a
						// node is on a crack
						if(warnings.Issue(CrackHeader::warnNodeOnCrack,11)==GAVE_WARNING) Describe();
						
						// tell field [0] it has crack from field 1 (but info currently not used)
						if(cvf[0]->location(FIRST_CRACK)==NO_CRACK)
							cvf[0]->SetLocationAndCrack(cfld[0].loc,cfld[0].crackNum,FIRST_CRACK);
					}
				}
						
				// Here means found a new crack, which hopefully will be appropriate for cvf[2]
				// Here means cvf[1]->crackNumber(FIRST_CRACK)!=cfld[0].crackNum
				else
				{	// create [2] is possible, otherwise see if same crack or not
					if(!CrackVelocityField::ActiveField(cvf[2]))
					{	if(cvf[2]==NULL)
						{	cvf[2]=CrackVelocityField::CreateCrackVelocityField(cfld[0].loc,cfld[0].crackNum);
							if(cvf[2]==NULL) throw CommonException("Memory error allocating crack velocity field 2.",
																   "NodalPoint::AddMomentumTask1");
						}
						else
							cvf[2]->SetLocationAndCrack(cfld[0].loc,cfld[0].crackNum,FIRST_CRACK);
						vfld=2;
					}
					else if(cvf[2]->crackNumber(FIRST_CRACK)==cfld[0].crackNum)
					{	if(cvf[2]->location(FIRST_CRACK)==cfld[0].loc)
						{	// found another point for [2]
							vfld=2;
						}
						else
						{	// it can only be field 0
							vfld=0;
						
							// Here means both above and below crack for field [2], which can only happen is a
							// node is on a crack
							if(warnings.Issue(CrackHeader::warnNodeOnCrack,12)==GAVE_WARNING) Describe();
							
							// tell field [0] it has crack from field 2 (but info currently not used)
							if(cvf[0]->location(SECOND_CRACK)==NO_CRACK)
								cvf[0]->SetLocationAndCrack(cfld[0].loc,cfld[0].crackNum,SECOND_CRACK);
						}
					}
					else
					{	// found a third crack at this node, try to use [0] and issue warning
						vfld=0;
						warnings.Issue(CrackHeader::warnThreeCracks,11);
					}
				}
				
				// add the crack normals to selected field
				cvf[vfld]->AddNormals(&cfld[0].norm,FIRST_CRACK);
				
				break;
			default:
				break;
		}
	}
	
	// two fields - currently only allows one field and thus cannot handle nodes on interacting cracks
	// two fields are always put into cvf[3] and it is only field with alternate crack information
	else
	{	if(!CrackVelocityField::ActiveField(cvf[3]))
		{	// store in any order
			if(cvf[3]==NULL)
			{	cvf[3]=CrackVelocityField::CreateCrackVelocityField(cfld[0].loc,cfld[0].crackNum);
				if(cvf[3]==NULL) throw CommonException("Memory error allocating crack velocity field 3.",
														  "NodalPoint::AddMomentumTask1");
			}
			else
				cvf[3]->SetLocationAndCrack(cfld[0].loc,cfld[0].crackNum,FIRST_CRACK);
			cvf[3]->SetLocationAndCrack(cfld[1].loc,cfld[1].crackNum,SECOND_CRACK);
		}
		else
		{	// make sure same two cracks
			int c1=-1,c2=-1;
			if(cfld[0].crackNum==cvf[3]->crackNumber(FIRST_CRACK) && cfld[1].crackNum==cvf[3]->crackNumber(SECOND_CRACK))
			{	c1=0;
				c2=1;
			}
			else if(cfld[0].crackNum==cvf[3]->crackNumber(SECOND_CRACK) && cfld[1].crackNum==cvf[3]->crackNumber(FIRST_CRACK))
			{	c1=1;
				c2=0;
			}
			if(c1>=0)
			{	if(cfld[c1].loc==cvf[3]->location(FIRST_CRACK))
				{	// first crack is correct orientation
					if(cfld[c2].loc!=cvf[3]->location(SECOND_CRACK))
					{	// second crack is wrong or node may be on second crack, might be better to use [1]
						if(warnings.Issue(CrackHeader::warnNodeOnCrack,21)==GAVE_WARNING) Describe();
					}
				}
				else if(cfld[c2].loc==cvf[3]->location(SECOND_CRACK))
				{	// second crack correct, but first is wrong or node may be on second crack, might be better to use [2]
					if(warnings.Issue(CrackHeader::warnNodeOnCrack,22)==GAVE_WARNING) Describe();
				}
				else
				{	// both cracks wrong, not sure what to do
					if(warnings.Issue(CrackHeader::warnNodeOnCrack,23)==GAVE_WARNING) Describe();
				}
			}
			else
			{	// at least a third crack found for this node, not sure what to do
				warnings.Issue(CrackHeader::warnThreeCracks,21);
			}
		}
		vfld=3;
		
		// add the crack normals to selected field
		cvf[vfld]->AddNormals(&cfld[0].norm,FIRST_CRACK);
		cvf[vfld]->AddNormals(&cfld[1].norm,SECOND_CRACK);
	}
	
    // add momemtum to selected field (in g mm/sec)
	Vector wtvel;
	cvf[vfld]->AddMomentumTask1(matfld,CopyScaleVector(&wtvel,vel,wt),vel);
	
	// return crack velocity field that was used
	return vfld;
}

// Add mass for selected field
void NodalPoint::AddVolumeGradient(short vfld,int matfld,MPMBase *mptr,double dNdx,double dNdy,double dNdz)
{	if(fmobj->multiMaterialMode)
		cvf[vfld]->AddVolumeGradient(matfld,mptr,dNdx,dNdy,dNdz);
}

// Calculate total mass. Calculations might need to exclude nodes whose
// mass is too small in crack calculations
void NodalPoint::CalcTotalMassAndCount(void)
{	int i;
	mass=0.;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
			mass+=cvf[i]->GetTotalMassAndCount();
	}
}

// When has rigid particles, multimaterial mode, and cracks, sum all rigid particles on
//	this node and transfer copy to all crack velocity fields
// This is only called if COMBINE_RIGID_MATERIALS is defined
void NodalPoint::CombineRigidParticles(void)
{	int i,j;
	
	// find first material velocity field for rigid particles
	// a node can only have one rigid material, otherwise the contact algorithm fails
	// thus finding first is enough
	MatVelocityField *rmvf=NULL;
	int rigidFieldNum;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
		{	rmvf=((CrackVelocityFieldMulti *)cvf[i])->GetRigidMaterialField(&rigidFieldNum);
			if(rmvf!=NULL) break;
		}
	}
	
	// if none, then done
	if(rmvf==NULL) return;
	
	// sum other fields with this same rigid material into that materials field in cvf[i]
	// since 0 to i-1 do not have rigid material, only need to look from i+1 to end
	for(j=i+1;j<maxCrackFields;j++)
	{	if(CrackVelocityField::ActiveField(cvf[j]))
		{	// copy rigid material from field j to field i
			((CrackVelocityFieldMulti *)cvf[i])->CombineRigidFrom((CrackVelocityFieldMulti *)cvf[j],rigidFieldNum);
		}
	}
	
	// transfer the summed result to all other fields
	for(j=0;j<maxCrackFields;j++)
	{	if(j!=i && CrackVelocityField::ActiveField(cvf[j]))
		{	// copy rigid material from field i to field ji
			((CrackVelocityFieldMulti *)cvf[j])->CopyRigidFrom((CrackVelocityFieldMulti *)cvf[i],rigidFieldNum);
		}
	}
}

#pragma mark TASK 3 METHODS

// Add to internal force
void NodalPoint::AddFintTask3(short vfld,int matfld,Vector *f) { cvf[vfld]->AddFintTask3(matfld,f); }

// Add to internal force spread out over materials for same acceleration on each
void NodalPoint::AddFintSpreadTask3(short vfld,Vector f) { cvf[vfld]->AddFintSpreadTask3(&f); }

// Add to external force (g-mm/sec^2)
void NodalPoint::AddFextTask3(short vfld,int matfld,Vector *f) { cvf[vfld]->AddFextTask3(matfld,f); }

// Add to external force spread out over materials for same acceleration on each
void NodalPoint::AddFextSpreadTask3(short vfld,Vector f) { cvf[vfld]->AddFextSpreadTask3(&f); }

// Calculate total force at a node from current values (now m*a in g mm/sec^2)
void NodalPoint::CalcFtotTask3(double extDamping)
{	int i;
    for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->CalcFtotTask3(extDamping);
	}
}

#pragma mark TASK 4 METHODS

// update momenta for this MPM step
void NodalPoint::UpdateMomentaOnNode(double timestep)
{	// update momenta
	int i;
    for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->UpdateMomentaOnField(timestep);
    }
}

#pragma mark TASK 5 METHODS

// Increment velocity and acceleration for this material point using one velocity field
void NodalPoint::IncrementDelvaTask5(short vfld,int matfld,double fi,Vector *delv,Vector *dela)
{	cvf[vfld]->IncrementDelvaTask5(matfld,fi,delv,dela);
}

#pragma mark TASK 6 METHODS

// zero momentum at a node for new calculations
void NodalPoint::RezeroNodeTask6(double deltaTime)
{	int i;
    for(i=0;i<maxCrackFields;i++)
    {	if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->RezeroNodeTask6(deltaTime);
    }
}

// Add 2 momentum on second pass for selected field that must be there
void NodalPoint::AddMomentumTask6(short vfld,int matfld,double wt,Vector *vel)
{	cvf[vfld]->AddMomentumTask6(matfld,wt,vel);
}

#pragma mark TASK 7 CRACK CALCULATIONS FOR J AND K

// Initialize fields for grid extrapolations for strains, etc.
void NodalPoint::ZeroDisp(void)
{	
	int i;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveNonrigidField(cvf[i]))
			cvf[i]->CreateStrainField();
	}
	 
	// Five possibitilies (s), (s,a), (s,b), (a,b), (b,a) where
	//	s is same side of crack, a and b are above and below.
	// Field [1], if present, tells which is above or below
	// This calculation assumes only 1 crack
	if(!CrackVelocityField::ActiveNonrigidField(cvf[1]))
	{	// only field [0] so both are zero
		above=below=0;
	}
	else if(cvf[1]->location(FIRST_CRACK)==ABOVE_CRACK)
	{	below=0;
		above=1;
	}
	else
	{	below=1;
		above=0;
	}
}

// delete any strain fields that were created
void NodalPoint::DeleteDisp(void)
{
	int i;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveNonrigidField(cvf[i]))
			cvf[i]->DeleteStrainField();
	}
}

// Add to displacement gradient
void NodalPoint::AddUGradient(short vfld,double wt,double dudx,double dudy,double dvdx,double dvdy)
{	DispField *df=cvf[vfld]->df;
	df->du.x+=wt*dudx;
	df->du.y+=wt*dudy;
	df->dv.x+=wt*dvdx;
	df->dv.y+=wt*dvdy;
}

// Add to kinetic energy and strain energy
// wt includes density (g/cm^3)
// Must scale kinetic later by 0.5e-9 to get units N/mm^2 (here is twice g/cm^3 mm^2/sec^2)
// Must scale work later by 1.e-6 to get N/mm^2 (here is N/m^2)
void NodalPoint::AddEnergy(short vfld,double wt,double vx,double vy,double work)
{	DispField *df=cvf[vfld]->df;
	df->kinetic+=wt*(vx*vx+vy*vy);
	df->work+=wt*work;
}

// add to nodal stress (in-plane only)
// wt includes rho thus stress is N/m^2
// Must scale by 1e-6 to get N/mm^2
void NodalPoint::AddStress(short vfld,double wt,Tensor *stress)
{	DispField *df=cvf[vfld]->df;
	df->stress.xx+=wt*stress->xx;
	df->stress.yy+=wt*stress->yy;
	df->stress.xy+=wt*stress->xy;
}

// Finish extrapolated strain field terms
//	Divide by mass to get final results
//	Scale to compatible units
void NodalPoint::CalcStrainField(void)
{
	int j;
	double mnode;

	// do all strain fields
	for(j=0;j<maxCrackFields;j++)
	{	if(!CrackVelocityField::ActiveNonrigidField(cvf[j])) continue;
		DispField *df=cvf[j]->df;
		if(df==NULL) continue;
		
		mnode=1./cvf[j]->GetTotalMass();
		df->du.x*=mnode;			// no units
		df->du.y*=mnode;
		df->dv.x*=mnode;
		df->dv.y*=mnode;
		df->kinetic*=mnode*.5e-9;	// N/mm^2
		mnode*=1.e-6;
		df->work*=mnode;			// N/mm^2
		df->stress.xx*=mnode;		// N/mm^2
		df->stress.yy*=mnode;
		df->stress.xy*=mnode;
	}
}

// interpolate two nodes (near crack plane). This method is only called in J calculation for the
// phatom node placed at the crack plane.
// Symbolically gets fract*n2 + (1-fract)*n1
void NodalPoint::Interpolate(NodalPoint *n1,NodalPoint *n2,double fract,bool startTip)
{
	// Dereference above and below fields for each crack
	DispField *a1fld=n1->cvf[(int)n1->above]->df;
	DispField *b1fld=n1->cvf[(int)n1->below]->df;
	DispField *a2fld=n2->cvf[(int)n2->above]->df;
	DispField *b2fld=n2->cvf[(int)n2->below]->df;
	
	// need strain field in first crack velocity field and entire second
	// crack velocity field for this phatom node (it may be zero)
	cvf[0]->CreateStrainField();
	cvf[1]=CrackVelocityField::CreateCrackVelocityField(0,0);
	if(cvf[1]==NULL) throw CommonException("Memory error allocating crack velocity field 1.",
										   "NodalPoint::Interpolate");
	cvf[1]->CreateStrainField();
	above=0;
	below=1;

	// node 2 has only one field
	if(n2->above==n2->below)
	{	if(startTip)
		{	// node 2 field is entirely above the crack (in fld [0]) or no below the crack field
			b2fld=NULL;
			if(n1->above==n1->below)
			{	// node 1 field is entirely below the crack (in fld [0]) or no above the crack field
				a1fld=NULL;
			}
		}
		else
		{	// node 2 field is entirely below the crack (in fld [0]) or no above the crack field
			a2fld=NULL;
			if(n1->above==n1->below)
			{	// node 1 field is entirely above the crack (in fld [0]) or no below the crack field
				b1fld=NULL;
			}
		}
	}

	// node 1 has only one field
	else if(n1->above==n1->below)
	{	if(startTip)
		{	// node 1 field is entirely below the crack (in fld [0]) or no above the crack field
			a1fld=NULL;
			// node 2 is mixed because was not trapped above
		}
		else
		{	// node 1 field is entirely above the crack (in fld [0]) or no below the crack field
			b1fld=NULL;
			// node 2 is mixed because was not trapped above
		}
	}

	// average final settings
	AverageStrain(cvf[(int)above]->df,a1fld,a2fld,fract);
	AverageStrain(cvf[(int)below]->df,b1fld,b2fld,fract);
}

// interpolate between two fields and store in destination field
void NodalPoint::AverageStrain(DispField *dest,DispField *src1,DispField *src2,double fract)
{	if(src1!=NULL && src2!=NULL)
	{	dest->du.x=fract*src2->du.x + (1.-fract)*src1->du.x;
		dest->du.y=fract*src2->du.y + (1.-fract)*src1->du.y;
		dest->dv.x=fract*src2->dv.x + (1.-fract)*src1->dv.x;
		dest->dv.y=fract*src2->dv.y + (1.-fract)*src1->dv.y;
		dest->kinetic=fract*src2->kinetic + (1.-fract)*src1->kinetic;
		dest->work=fract*src2->work + (1.-fract)*src1->work;
		dest->stress.xx=fract*src2->stress.xx + (1.-fract)*src1->stress.xx;
		dest->stress.yy=fract*src2->stress.yy + (1.-fract)*src1->stress.yy;
		dest->stress.zz=fract*src2->stress.zz + (1.-fract)*src1->stress.zz;
		dest->stress.xy=fract*src2->stress.xy + (1.-fract)*src1->stress.xy;
	}
	else if(src1!=NULL)
	{	dest->du.x=(1.-fract)*src1->du.x;
		dest->du.y=(1.-fract)*src1->du.y;
		dest->dv.x=(1.-fract)*src1->dv.x;
		dest->dv.y=(1.-fract)*src1->dv.y;
		dest->kinetic=(1.-fract)*src1->kinetic;
		dest->work=(1.-fract)*src1->work;
		dest->stress.xx=(1.-fract)*src1->stress.xx;
		dest->stress.yy=(1.-fract)*src1->stress.yy;
		dest->stress.zz=(1.-fract)*src1->stress.zz;
		dest->stress.xy=(1.-fract)*src1->stress.xy;
	}
	else if(src2!=NULL)
	{	dest->du.x=fract*src2->du.x;
		dest->du.y=fract*src2->du.y;
		dest->dv.x=fract*src2->dv.x;
		dest->dv.y=fract*src2->dv.y;
		dest->kinetic=fract*src2->kinetic;
		dest->work=fract*src2->work;
		dest->stress.xx=fract*src2->stress.xx;
		dest->stress.yy=fract*src2->stress.yy;
		dest->stress.zz=fract*src2->stress.zz;
		dest->stress.xy=fract*src2->stress.xy;
	}
	else
	{	// outside the grid
		dest->du.x=0.;
		dest->du.y=0.;
		dest->dv.x=0.;
		dest->dv.y=0.;
		dest->kinetic=0.;
		dest->work=0.;
		dest->stress.xx=0.;
		dest->stress.yy=0.;
		dest->stress.zz=0.;
		dest->stress.xy=0.;
	}
}

#pragma mark TASK 8 METHODS

// Increment velocity for crack surface
short NodalPoint::IncrementDelvSideTask8(short side,int crackNumber,double fi,Vector *delv,double *surfaceMass,CrackSegment *seg)
{
	short vfld=-1;
	
	double x1=seg->surfx[side-1];
	double y1=seg->surfy[side-1];
	
	// check in field [1] if it is for crack crackNumber
	if(CrackVelocityField::ActiveNonrigidField(cvf[1]))
	{	if(cvf[1]->crackNumber(FIRST_CRACK)==crackNumber)
		{	if(side==cvf[1]->location(FIRST_CRACK))
			{	vfld=1;
				
				// maybe switch [1] to [3]
				if(CrackVelocityField::ActiveNonrigidField(cvf[3]))
				{	// if line crosses second crack in [3], switch to [3]
					int otherCrack=cvf[3]->OppositeCrackTo(crackNumber,side);
					if(otherCrack>0)
					{	if(SurfaceCrossesOneCrack(x1,y1,x,y,otherCrack)!=NO_CRACK)
							vfld=3;
					}
				}
				else if(SurfaceCrossesOtherCrack(x1,y1,x,y,crackNumber))
					vfld=3;
			}
			else
			{	vfld=0;
				
				// maybe switch [0] to [2]
				if(CrackVelocityField::ActiveNonrigidField(cvf[2]))
				{	// if line crosses the crack found in [2], then switch to [2]
					if(SurfaceCrossesOneCrack(x1,y1,x,y,cvf[2]->crackNumber(FIRST_CRACK))!=NO_CRACK)
					   vfld=2;
				}
				else if(SurfaceCrossesOtherCrack(x1,y1,x,y,crackNumber))
					vfld=2;
			}
		}
	}
	
	// if not found in [1], check in [2]
	if(vfld<0 && CrackVelocityField::ActiveNonrigidField(cvf[2]))
	{	if(cvf[2]->crackNumber(FIRST_CRACK)==crackNumber)
		{	if(side==cvf[2]->location(FIRST_CRACK))
			{	vfld=2;
		
				// maybe switch [2] to [3]
				if(CrackVelocityField::ActiveNonrigidField(cvf[3]))
				{	// if line crosses second crack in [3], switch to [3]
					int otherCrack=cvf[3]->OppositeCrackTo(crackNumber,side);
					if(otherCrack>0)
					{	if(SurfaceCrossesOneCrack(x1,y1,x,y,otherCrack)!=NO_CRACK)
							vfld=3;
					}
				}
				else if(SurfaceCrossesOtherCrack(x1,y1,x,y,crackNumber))
					vfld=3;
			}
			else
			{	vfld=0;
				
				// maybe switch [0] to [1]
				if(CrackVelocityField::ActiveNonrigidField(cvf[1]))
				{	// if line crosses the crack in [1], then switch to [1]
					if(SurfaceCrossesOneCrack(x1,y1,x,y,cvf[1]->crackNumber(FIRST_CRACK))!=NO_CRACK)
						vfld=1;
				}
				else if(SurfaceCrossesOtherCrack(x1,y1,x,y,crackNumber))
					vfld=1;
			}
		}
	}
	
	// if not found in [1] or [2], look in [3]
	if(vfld<0 && CrackVelocityField::ActiveNonrigidField(cvf[3]))
	{	// verify has correct field and retreive other crack
		int otherCrack=cvf[3]->OppositeCrackTo(crackNumber,side);
		if(otherCrack>0)
		{	// if crack found, then can use [3] is crosses the other crack in [3]
			if(SurfaceCrossesOneCrack(x1,y1,x,y,otherCrack)!=NO_CRACK)
				vfld=3;
		}
		else
		{	// when not found, try to see if cross path other crack to [1] or [2]
			otherCrack=cvf[3]->OppositeCrackTo(crackNumber,ABOVE_CRACK+BELOW_CRACK-side);
			if(otherCrack>0)
			{	if(CrackVelocityField::ActiveNonrigidField(cvf[1]))
				{	if(otherCrack==cvf[1]->crackNumber(FIRST_CRACK))
					{	if(SurfaceCrossesOneCrack(x1,y1,x,y,otherCrack)!=NO_CRACK)
							vfld=1;
						else
							vfld=0;
					}
				}
				if(vfld<0 && CrackVelocityField::ActiveNonrigidField(cvf[2]))
				{	if(otherCrack==cvf[2]->crackNumber(FIRST_CRACK))
					{	if(SurfaceCrossesOneCrack(x1,y1,x,y,otherCrack)!=NO_CRACK)
							vfld=2;
						else
							vfld=0;
					}
				}
				// if fails, will try [0] below
			}
		}
	}
	
	// if still not found, see if [0] can be used
	if(vfld<0 && CrackVelocityField::ActiveNonrigidField(cvf[0]))
	{	Vector moved=seg->SlightlyMoved(side);
		CrackField cfld[2];
		SurfaceCrossesCracks(moved.x,moved.y,x,y,cfld);
        //cout << "#Active on node " << num << ": " << CrackVelocityField::ActiveNonrigidField(cvf[0]) << "," <<
        //                       CrackVelocityField::ActiveNonrigidField(cvf[1]) << endl;
        //cout << "#Side " << side << " (" << x1 << "," << y1 << ") to (" << moved.x << "," << moved.y << ")" << endl;
		if(cfld[0].loc==NO_CRACK)
			vfld=0;
		else if(cfld[1].loc==NO_CRACK)
		{	// only one crack was found - does it match [1] or [2]
			if(CrackVelocityField::ActiveNonrigidField(cvf[1]))
			{	if(cfld[0].crackNum==cvf[1]->crackNumber(FIRST_CRACK))
				{	if(cfld[0].loc==cvf[1]->location(FIRST_CRACK))
						vfld=1;
					else
						vfld=0;				// if surface particle was on the crack
				}
			}
			if(vfld<0 && CrackVelocityField::ActiveNonrigidField(cvf[2]))
			{	if(cfld[0].crackNum==cvf[2]->crackNumber(FIRST_CRACK))
				{	if(cfld[0].loc==cvf[2]->location(FIRST_CRACK))
						vfld=2;
					else
						vfld=0;				// if surface particle was on the crack
				}
			}
		}
		else if(CrackVelocityField::ActiveNonrigidField(cvf[3]))
		{	// found two cracks, but only use if same two cracks that are in [3]
			if(cfld[0].crackNum==cvf[3]->crackNumber(FIRST_CRACK) && cfld[1].crackNum==cvf[3]->crackNumber(SECOND_CRACK))
			{	if(cfld[0].loc==cvf[3]->location(FIRST_CRACK) && cfld[1].loc==cvf[3]->location(SECOND_CRACK))
					vfld=3;
			}
			else if(cfld[1].crackNum==cvf[3]->crackNumber(FIRST_CRACK) && cfld[0].crackNum==cvf[3]->crackNumber(SECOND_CRACK))
			{	if(cfld[1].loc==cvf[3]->location(FIRST_CRACK) && cfld[0].loc==cvf[3]->location(SECOND_CRACK))
					vfld=3;
			}
		}
	}
	
	// Exit if no field
	if(vfld<0) return FALSE;
	if(!CrackVelocityField::ActiveNonrigidField(cvf[vfld])) return FALSE;
	
	// increment the velocity if enough mass
	double fieldMass=mass;
	if(cvf[vfld]->IncrementDelvTask8(fi,delv,&fieldMass))
	{	*surfaceMass+=fi*fieldMass;
		return TRUE;
	}
	return FALSE;
}

// Determine if line from crack surface particle (1) to node (2) crosses 
// one or more cracks
void NodalPoint::SurfaceCrossesCracks(double x1,double y1,double x2,double y2,CrackField *cfld)
{
	CrackHeader *nextCrack=firstCrack;
	int cfound=0;
    short vfld;
	Vector norm;
	
	cfld[0].loc=NO_CRACK;			// NO_CRACK, ABOVE_CRACK, or BELOW_CRACK
	cfld[1].loc=NO_CRACK;
	
	while(nextCrack!=NULL)
	{	vfld=nextCrack->CrackCross(x1,y1,x2,y2,&norm);
		if(vfld!=NO_CRACK)
		{	cfld[cfound].loc=vfld;
			cfld[cfound].crackNum=nextCrack->GetNumber();
			cfound++;
			if(cfound>1) break;			// stop if found two, if there are more then two, physics may be off
		}
		nextCrack=(CrackHeader *)nextCrack->GetNextObject();
	}
}

// Determine if line from crack surface particle (1) to node (2) crosses specific crack
// one or more cracks
int NodalPoint::SurfaceCrossesOneCrack(double x1,double y1,double x2,double y2,int cnum)
{
	CrackHeader *nextCrack=firstCrack;
	Vector norm;
	
	while(nextCrack!=NULL)
	{	if(nextCrack->GetNumber()==cnum)
			return nextCrack->CrackCross(x1,y1,x2,y2,&norm);
		nextCrack=(CrackHeader *)nextCrack->GetNextObject();
	}
	return NO_CRACK;
}

// Determine if line from crack surface particle (1) to node (2) crosses any crack other than the
// one provided
bool NodalPoint::SurfaceCrossesOtherCrack(double x1,double y1,double x2,double y2,int cnum)
{
	CrackHeader *nextCrack=firstCrack;
	Vector norm;
	
	while(nextCrack!=NULL)
	{	if(nextCrack->GetNumber()!=cnum)
		{	if(nextCrack->CrackCross(x1,y1,x2,y2,&norm)!=NO_CRACK)
				return TRUE;
		}
		nextCrack=(CrackHeader *)nextCrack->GetNextObject();
	}
	return FALSE;
}

// Calculate CM velocity at a node and store in nv[0] velocity
//   (only used when contact.GetMoveOnlySurfaces() is FALSE, i.e., when crack particles
//		move in CM velocity field,  and then only when crack plane particles are about to move)
void NodalPoint::CalcCMVelocityTask8(void)
{
	int i;
	Vector nodePk;
	ZeroVector(&nodePk);
	int totalParticles=0;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveNonrigidField(cvf[i]))
			totalParticles+=cvf[i]->CollectMomentaTask8(&nodePk);
	}
	
	// store in field zero
	cvf[0]->SetCMVelocityTask8(ScaleVector(&nodePk, 1./mass),totalParticles);
}

// Get velocity for center of mass
//   (assumes recently called CalcCMVelocityTask8() for this node)
//   (only used when contact.GetMoveOnlySurfaces() is FALSE and crack plane
//       particles are moving)
bool NodalPoint::GetCMVelocityTask8(Vector *vk)
{	return cvf[0]->GetCMVelocityTask8(vk);
}


#pragma mark INCREMENTERS

// Add displacements to selected field, but in task 6, add rigid particle velocity
// to all crack fields (will have vfld<0)
void NodalPoint::AddDisplacement(short vfld,int matfld,double wt,Vector *pdisp)
{	cvf[vfld]->AddDisplacement(matfld,wt,pdisp);
}

// Add volume to selected field
void NodalPoint::AddVolume(short vfld,int matfld,double wtVol)
{	cvf[vfld]->AddVolume(matfld,wtVol);
}

#pragma mark VELOCITY FIELDS

// Calculate velocity at a node from current momentum and mass matrix
void NodalPoint::CalcVelocityForStrainUpdate(void)
{	// get velocity for all crack and material velocity fields
	int i;
    for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->CalcVelocityForStrainUpdate();
    }
}

// Get velocity for selected field
Vector NodalPoint::GetVelocity(short vfld,int matfld)
{	return cvf[vfld]->GetVelocity(matfld);
}

// Get velocity for selected field
Vector NodalPoint::GetContactForce(short vfld,int matfld)
{	return cvf[vfld]->GetContactForce(matfld);
}

// Some all forces from rigid material velocity fields
Vector NodalPoint::GetTotalContactForce(bool clearForces)
{	int i;
	Vector fcontact;
	ZeroVector(&fcontact);
    for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->SumAndClearRigidContactForces(&fcontact,clearForces);
	}
	double scale=-1.e-6/timestep;
	fcontact.x*=scale;
	fcontact.y*=scale;
	fcontact.z*=scale;
	return fcontact;
}

#pragma mark MATERIAL CONTACT

// Called in multimaterial mode to check contact at nodes with multiple materials
void NodalPoint::MaterialContactOnNode(bool postUpdate,double deltime)
{
	// check each crack velocity field on this node
	for(int i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->MaterialContact(num,i,postUpdate,deltime);
	}
}

// retrieve -2*scale*(mass gradient) for material matfld in velocity field vfld
void NodalPoint::GetVolumeGradient(short vfld,int matfld,Vector *grad,double scale)
{	cvf[vfld]->GetVolumeGradient(matfld,this,grad,scale);
}

// This node is known to have imperfect interface with forces in cvf[vfld] from material mati
// to matipaired (second material with the max volume)
void NodalPoint::MaterialInterfaceForce(MaterialInterfaceNode *mmnode)
{	
    // recall interface response from material interface node
    Vector fImpInt;
    double energy = mmnode->GetInterfaceTraction(&fImpInt);
    
    // add total force (in g mm/sec^2) to material field
    int vfld,mati,matipaired;
    mmnode->GetFieldInfo(&vfld, &mati, &matipaired);
    cvf[vfld]->AddFintTask3(mati,&fImpInt);
    
    // add negative force to paired material (in a pair)
    if(matipaired>=0)
    {   ScaleVector(&fImpInt, -1.);
        cvf[vfld]->AddFintTask3(matipaired,&fImpInt);
    }
    
    // add interface energy in units g-mm^2/sec^2 (multiply by 1e-9 to get J - kg-m^2/sec^2)
    interfaceEnergy+=energy;
}

#pragma mark CRACK SURFACE CONTACT

// Look for crack contact and adjust accordingly.
void NodalPoint::CrackContact(int makeCopy,bool postUpdate,double deltime)
{	// Nothing to do if not near a crack contact surface: Possible fields are
	//  1. Those with no contacts: [0], [1], [3], [0]&[3], [1]&[2]
	//  2. Those with contacts: [0]&[1], [1]&[3], [0]&[1]&[2], [0]&[1]&[3], [1]&[2]&[3], and [0]&[1]&[2]&[3]
	//  3. Never occurs [2], [0]&[2], [2]&[3], [0]&[2]&[3]
	
	// exit on no contact
	bool has1=CrackVelocityField::ActiveNonrigidField(cvf[1]);
	bool has2=CrackVelocityField::ActiveNonrigidField(cvf[2]);
	if(!has1 && !has2) return;	// True for [0], [3], and [0]&[3]
	bool has0=CrackVelocityField::ActiveNonrigidField(cvf[0]);
	bool has3=CrackVelocityField::ActiveNonrigidField(cvf[3]);
	if(!has0 && !has3) return;	// True for [1] and [1]&[2]
	
	// store references to this node for future use
	if(makeCopy)
	{	CrackNode::currentNode=new CrackNode(this);
		if(CrackNode::currentNode==NULL) throw CommonException("Memory error allocating storage for a crack node.",
															   "NodalPoint::CrackContact");
	}
	
	// between [0] and [1]
	int cnum,cabove;
	if(has0 && has1)
	{	cnum=cvf[1]->crackNumber(FIRST_CRACK);
		if(contact.HasContact(cnum))
		{	cabove=(cvf[1]->location(FIRST_CRACK)==ABOVE_CRACK) ? 1 : 0;
			AdjustContact(cabove,1-cabove,&(cvf[1]->norm[FIRST_CRACK]),cnum,postUpdate,deltime);
		}
	}
	
	// between [0] & [2]
	if(has0 && has2)
	{	cnum=cvf[2]->crackNumber(FIRST_CRACK);
		if(contact.HasContact(cnum))
		{	cabove=(cvf[2]->location(FIRST_CRACK)==ABOVE_CRACK) ? 2 : 0;
			AdjustContact(cabove,2-cabove,&(cvf[2]->norm[FIRST_CRACK]),cnum,postUpdate,deltime);
		}
	}
	
	// between [1] & [3]
	if(has1 && has3) CrackContactThree(1,postUpdate,deltime);

	// between [2] & [3]
	if(has2 && has3) CrackContactThree(2,postUpdate,deltime);
	
}

// Contact between field [single] and field [3] when both fields are present
// Possibilities
//	1. [single] and [3] opposite sides of crack for [single], do contact that crack
//	2. [single] and [3] same side of crack for single, then do contact across other
//			crack form [3]
//	3. Neither crack in [3] is same as crack for [single], then no contact
void NodalPoint::CrackContactThree(int single,bool postUpdate,double deltime)
{
	int cnum=cvf[single]->crackNumber(FIRST_CRACK);
	int cloc=cvf[single]->location(FIRST_CRACK);
	
	int matchCrack;
	if(cnum==cvf[3]->crackNumber(FIRST_CRACK))
		matchCrack=FIRST_CRACK;
	else if(cnum==cvf[3]->crackNumber(SECOND_CRACK))
		matchCrack=SECOND_CRACK;
	else
		return;

	int cabove;
	Vector cnorm;
	if(cloc!=cvf[3]->location(matchCrack))
	{	// [single] and [3] are on opposite sides of cnum
		cabove = (cloc==ABOVE_CRACK) ? single : 3 ;
		cnorm=cvf[single]->norm[FIRST_CRACK];
	}
	else
	{	// [single] and [3] same side as cnum, do contact across the other crack
		matchCrack = (matchCrack==FIRST_CRACK) ? SECOND_CRACK : FIRST_CRACK ;
		cnum=cvf[3]->crackNumber(matchCrack);
		cabove = (cvf[3]->crackNumber(matchCrack)) ? 3 : single ;
		ZeroVector(&cnorm);
	}
	
	// do contact if this crack has contact
	if(contact.HasContact(cnum))
	{	AddVector(&cnorm,&(cvf[3]->norm[matchCrack]));			// include [3] normal
		AdjustContact(cabove,single+3-cabove,&cnorm,cnum,postUpdate,deltime);
	}
}

// Look for crack contact and adjust accordingly - a for field above and b for field below and both
// fields must be verified as present (1 or more points)
void NodalPoint::AdjustContact(short a,short b,Vector *norm,int crackNumber,bool postUpdate,double deltime)
{
#ifdef _BC_CRACK_SIDE_ONLY_
	Vector delPa,delPb;
	if(!contact.GetDeltaMomentum(this,&delPa,&delPb,cvf[a],cvf[b],norm,crackNumber,postUpdate,deltime))
		return;
	
	// on post update contact, do not change nodes with boundary conditions
	if(postUpdate && fixedDirection)
	{	if(fixedDirection&X_DIRECTION)
		{	delPa.x=0.;
			delPb.x=0.;
		}
		if(fixedDirection&Y_DIRECTION)
		{	delPa.y=0.;
			delPb.y=0.;
		}
		if(fixedDirection&Z_DIRECTION)
		{	delPa.z=0.;
			delPb.z=0.;
		}
	}
	
    // change momenta
	cvf[a]->ChangeMomentum(&delPa,postUpdate,deltime);
    cvf[b]->ChangeMomentum(&delPb,postUpdate,deltime);
#else
    Vector delP;
	if(!contact.GetDeltaMomentum(this,&delP,cvf[a],cvf[b],norm,crackNumber,postUpdate,deltime))
		return;
	
	// on post update contact, do not change nodes with boundary conditions
	if(postUpdate && fixedDirection)
	{	if(fixedDirection&X_DIRECTION) delP.x=0.;
		if(fixedDirection&Y_DIRECTION) delP.y=0.;
		if(fixedDirection&Z_DIRECTION) delP.z=0.;
	}
	
    // change momenta
	cvf[a]->ChangeMomentum(&delP,postUpdate,deltime);
	Vector delPb;
    cvf[b]->ChangeMomentum(CopyScaleVector(&delPb,&delP,-1.),postUpdate,deltime);
#endif
}

// Look for crack contact and adjust accordingly
void NodalPoint::CrackInterfaceForce(void)
{	// Nothing to do if not near a crack contact surface: Possible fields are
	//  1. Those with no contacts: [0], [1], [3], [0]&[3], [1]&[2]
	//  2. Those with contacts: [0]&[1], [1]&[3], [0]&[1]&[2], [0]&[1]&[3], [1]&[2]&[3], and [0]&[1]&[2]&[3]
	//  3. Never occurs [2], [0]&[2], [2]&[3], [0]&[2]&[3]
	
	// skip those with no contact
	bool has1=CrackVelocityField::ActiveNonrigidField(cvf[1]);
	bool has2=CrackVelocityField::ActiveNonrigidField(cvf[2]);
	if(!has1 && !has2) return;		// True for [0], [3], and [0]&[3]
	bool has0=CrackVelocityField::ActiveNonrigidField(cvf[0]);
	bool has3=CrackVelocityField::ActiveNonrigidField(cvf[3]);
	if(!has0 && !has3) return;	// True for [1] and [1]&[2]
    
	
	// between [0] and [1] across first crack
	int cnum,cabove;
	if(has0 && has1)
	{	cnum=cvf[1]->crackNumber(FIRST_CRACK);
		if(contact.IsImperfect(cnum))
		{	cabove=(cvf[1]->location(FIRST_CRACK)==ABOVE_CRACK) ? 1 : 0;
			AddInterfaceForce(cabove,1-cabove,&(cvf[1]->norm[FIRST_CRACK]),cnum);
		}
	}
			
	// between [0] and [2] across second crack
	if(has0 && has2)
	{	cnum=cvf[2]->crackNumber(FIRST_CRACK);
		if(contact.IsImperfect(cnum))
		{	cabove=(cvf[2]->location(FIRST_CRACK)==ABOVE_CRACK) ? 2 : 0;
			AddInterfaceForce(cabove,2-cabove,&(cvf[2]->norm[FIRST_CRACK]),cnum);
		}
	}
		
	
	// between [1] & [3]
	if(has1 && has3) InterfaceForceThree(1);
	
	// between [2] & [3]
	if(has2 && has3) InterfaceForceThree(2);
}

// Interface force between field [single] and field [3] when both fields are present
// Possibilities
//	1. [single] and [3] opposite sides of crack for [single], do force that crack
//	2. [single] and [3] same side of crack for single, then do force across other
//			crack from [3]
//	3. Neither crack in [3] is same as crack for [single], then no force
void NodalPoint::InterfaceForceThree(int single)
{
	int cnum=cvf[single]->crackNumber(FIRST_CRACK);
	int cloc=cvf[single]->location(FIRST_CRACK);
	
	int matchCrack;
	if(cnum==cvf[3]->crackNumber(FIRST_CRACK))
		matchCrack=FIRST_CRACK;
	else if(cnum==cvf[3]->crackNumber(SECOND_CRACK))
		matchCrack=SECOND_CRACK;
	else
		return;
	
	int cabove;
	Vector cnorm;
	if(cloc!=cvf[3]->location(matchCrack))
	{	// [single] and [3] are on opposite sides of cnum
		cabove = (cloc==ABOVE_CRACK) ? single : 3 ;
		cnorm=cvf[single]->norm[FIRST_CRACK];
	}
	else
	{	// [single] and [3] same side as cnum, do contact across the other crack
		matchCrack = (matchCrack==FIRST_CRACK) ? SECOND_CRACK : FIRST_CRACK ;
		cnum=cvf[3]->crackNumber(matchCrack);
		cabove = (cvf[3]->crackNumber(matchCrack)) ? 3 : single ;
		ZeroVector(&cnorm);
	}
	
	// do contact if this crack has contact
	if(contact.IsImperfect(cnum))
	{	AddVector(&cnorm,&(cvf[3]->norm[matchCrack]));			// include [3] normal
		AddInterfaceForce(cabove,single+3-cabove,&cnorm,cnum);
	}
}

// Look for cracks as imperfect interfaces and adjust accordingly - a for field above and b for field below
// fields must be verified as present (1 or more points)
void NodalPoint::AddInterfaceForce(short a,short b,Vector *norm,int crackNumber)
{	Vector fImpInt;
	
    // Use contact laws to change momenta - returns TRUE or FALSE if adjustment was made
	double rawEnergy;
	if(!contact.GetInterfaceForce(this,&fImpInt,cvf[a],cvf[b],norm,crackNumber,&rawEnergy))
		return;
	
    // Angled path correction method 1: distance to elipse through cell corners
    double dxnx = mpmgrid.gridx*norm->x, dyny = mpmgrid.gridy*norm->y ;
    double dist = sqrt((dxnx*dxnx+dyny*dyny)/(norm->x*norm->x+norm->y*norm->y));
	
	// Angled path correction method 2 (in imperfect interface by cracks paper):
    //   Find perpendicular distance which gets smaller as interface tilts
	//   thus the surface area increases
    //double dxnx = fabs(mpmgrid.gridx*norm->x), dyny = fabs(mpmgrid.gridy*norm->y) ;
	//double dist = fmax(dxnx,dyny)/sqrt(norm->x*norm->x+norm->y*norm->y);
    
	// Area correction method 1 (new): sqrt(2*vmin/vtot)*vtot/dist = sqrt(2*vmin*vtot)/dist
	double vola=cvf[a]->GetVolumeNonrigid(),volb=cvf[b]->GetVolumeNonrigid(),voltot=vola+volb;
	double surfaceArea=sqrt(2.0*fmin(vola,volb)*voltot)/dist;
	
	// Area correction method 2 (in imperfect interface by cracks paper): (2*vmin/vtot)*vtot/dist = 2*vmin/dist
	//double surfaceArea=2.0*fmin(cvf[a]->UnscaledVolumeNonrigid(),cvf[b]->UnscaledVolumeNonrigid())/dist;
    
    // If axisymmetric, multiply by radial position (vola, volb above were areas)
    if(fmobj->IsAxisymmetric()) surfaceArea *= x;
	
	// add total force (in g mm/sec^2)
	AddFintSpreadTask3(a,MakeVector(fImpInt.x*surfaceArea,fImpInt.y*surfaceArea,0.));
	AddFintSpreadTask3(b,MakeVector(-fImpInt.x*surfaceArea,-fImpInt.y*surfaceArea,0.));
	
	// add interface energy in units g-mm^2/sec^2 (multiply by 1e-9 to get J - kg-m^2/sec^2)
	interfaceEnergy+=rawEnergy*surfaceArea;
}

#pragma mark ACCESSORS

// number of particles for this node
int NodalPoint::NumberParticles(void)
{	int totalParticles=0;
	int i;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
			totalParticles+=cvf[i]->GetNumberPoints();
	}
	return totalParticles;
}

// number of particles for this node
int NodalPoint::NumberNonrigidParticles(void)
{	int totalParticles=0;
	int i;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
			totalParticles+=cvf[i]->GetNumberPointsNonrigid();
	}
	return totalParticles;
}

// describe velocity field
void NodalPoint::Describe(void)
{	cout << "# node=" << num << " pt=(" << x << "," << y << "," << z << ") mass=" << mass << endl;
	
	cout << "#  active crack velocity fields:" << endl;
	int i;
	int totalParticles=0,numFields=0;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(cvf[i]))
		{	cout << "#  " << i << ". ";
			cvf[i]->Describe();
			totalParticles+=cvf[i]->GetNumberPoints();
			numFields++;
		}
	}
}
	
#pragma mark BOUNDARY CONDITION METHODS

// set X velocity and momentum to zero
void NodalPoint::SetXMomVel(void)
{
#ifdef _BC_CRACK_SIDE_ONLY_
	// just set if on same side of crack
	cvf[0]->SetXMomVel();
#else
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->SetXMomVel();
	}
#endif
}

// Set Y velocity and momentum to zero
void NodalPoint::SetYMomVel(void)
{	
#ifdef _BC_CRACK_SIDE_ONLY_
	// just set if on same side of crack
	cvf[0]->SetYMomVel();
#else
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->SetYMomVel();
	}
#endif
}

// Set Z velocity and momentum to zero
void NodalPoint::SetZMomVel(void)
{	
#ifdef _BC_CRACK_SIDE_ONLY_
	// just set if on same side of crack
	cvf[0]->SetZMomVel();
#else
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->SetZMomVel();
	}
#endif
}

// Set velocity and momentum to zero in skewed direction (angle in radians)
// ccw from positive x axis
void NodalPoint::SetSkewMomVel(double angle)
{	
	// need to do skew condition at the velocity field level
#ifdef _BC_CRACK_SIDE_ONLY_
	// just set if on same side of crack
	cvf[0]->SetSkewMomVel(angle);
#else
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->SetSkewMomVel(angle);
	}
#endif
}

// Add X velocity and momentum at a node (assumes mass already set)
void NodalPoint::AddXMomVel(double velx)
{	
#ifdef _BC_CRACK_SIDE_ONLY_
	// just set if on same side of crack
	cvf[0]->AddXMomVelvelx);
#else
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->AddXMomVel(velx);
	}
#endif
}

// Add Y velocity and momentum at a node (assumes mass set)
void NodalPoint::AddYMomVel(double vely)
{
#ifdef _BC_CRACK_SIDE_ONLY_
	// just set if on same side of crack
	cvf[0]->AddYMomVel(vely);
#else
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->AddYMomVel(vely);
	}
#endif
}

// Add Z velocity and momentum at a node (assumes mass set)
void NodalPoint::AddZMomVel(double velz)
{	
#ifdef _BC_CRACK_SIDE_ONLY_
	// just set if on same side of crack
	cvf[0]->AddZMomVel(velz);
#else
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->AddZMomVel(velz);
	}
#endif
}

// Add velocity and momentum in skewed direction (angle in radians)
void NodalPoint::AddSkewMomVel(double vel,double angle)
{	
	// need to do skew condition at the velocity field level
#ifdef _BC_CRACK_SIDE_ONLY_
	// just set if on same side of crack
	cvf[0]->AddSkewMomVel(vel,angle);
#else
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->AddSkewMomVel(vel,angle);
	}
#endif
}

// set X force to -p(interpolated)/time such that updated momentum
//    of pk.x + deltime*ftot.x will be zero
void NodalPoint::SetXFtot(double deltime)
{	
#ifdef _BC_CRACK_SIDE_ONLY_
	// just on same side of the crack
	cvf[0]->SetXFtot(deltime);
#else
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->SetXFtot(deltime);
	}
#endif
}

// set Y force to -p(interpolated)/time such that updated momentum
//    of pk.y + deltime*ftot.y will be zero
void NodalPoint::SetYFtot(double deltime)
{	
#ifdef _BC_CRACK_SIDE_ONLY_
	// just on same side of the crack
	cvf[0]->SetYFtot(deltime);
#else
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->SetYFtot(deltime);
	}
#endif
}

// set Z force to -p(interpolated)/time such that updated momentum
//    of pk.z + deltime*ftot.z will be zero
void NodalPoint::SetZFtot(double deltime)
{	
#ifdef _BC_CRACK_SIDE_ONLY_
	// just on same side of the crack
	cvf[0]->SetZFtot(deltime);
#else
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->SetZFtot(deltime);
	}
#endif
}

// set skew force to -p(interpolated)/time such that updated momentum
//    of nv[j]->pk.x + deltime*nv[j]->ftot.x
//    of nv[j]->pk.y + deltime*nv[j]->ftot.y
// will be zero in skew direction
void NodalPoint::SetSkewFtot(double deltime,double angle)
{	
#ifdef _BC_CRACK_SIDE_ONLY_
	// just on same side of the crack
	cvf[0]->SetSkewFtot(deltime,angle);
#else
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->SetSkewFtot(deltime,angle);
	}
#endif
}

// set x force such that updated momentum will be mass*velocity
void NodalPoint::AddXFtot(double deltime,double velx)
{	
#ifdef _BC_CRACK_SIDE_ONLY_
	// just on same side of the crack
	cvf[0]->AddXFtot(deltime,velx);
#else
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->AddXFtot(deltime,velx);
	}
#endif
}

// set y force suce that updated momentum will be mass*velocity
void NodalPoint::AddYFtot(double deltime,double vely)
{	
#ifdef _BC_CRACK_SIDE_ONLY_
	// just on same side of the crack
	cvf[0]->AddYFtot(deltime,vely);
#else
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->AddYFtot(deltime,vely);
	}
#endif
}

// set z force such that updated momentum will be mass*velocity
void NodalPoint::AddZFtot(double deltime,double velz)
{	
#ifdef _BC_CRACK_SIDE_ONLY_
	// just on same side of the crack
	cvf[0]->AddZFtot(deltime,velz);
#else
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->AddZFtot(deltime,velz);
	}
#endif
}

// set Y force suce that updated momentum will be mass*velocity
void NodalPoint::AddSkewFtot(double deltime,double vel,double angle)
{	
#ifdef _BC_CRACK_SIDE_ONLY_
	// just on same side of the crack
	cvf[0]->AddSkewFtot(deltime,vel,angle);
#else
	int i;
	for(i=0;i<maxCrackFields;i++)
	{   if(CrackVelocityField::ActiveField(cvf[i]))
			cvf[i]->AddSkewFtot(deltime,vel,angle);
	}
#endif
}

// Mark a direction as fixed by velocity BC
// Assume 1 means x, 2 means y, 4 means z, x+y (or 3) is skewed condition
void NodalPoint::SetFixedDirection(int dir)
{	fixedDirection|=dir;
}

// Unmark a direction as fixed by velocity BC
// Assume 1 means x, 2 means y, 4 means z, x+y (or 3) is skewed condition
void NodalPoint::UnsetFixedDirection(int dir)
{	fixedDirection^=dir;
}

#pragma mark CLASS METHODS

// zero all velocity fields at start of time step
void NodalPoint::PreliminaryCalcs(void)
{	int i;
    for(i=1;i<=nnodes;i++)
        nd[i]->PrepareForFields();
}

// adjust momenta at overlaping material velocity fields
// multimaterial mode only
void NodalPoint::MaterialContact(bool multiMaterials,bool postUpdate,double deltime)
{	if(!multiMaterials) return;
	
	// implement material contact here
	int i;
	for(i=1;i<=nnodes;i++)
	{	// Each node calls cvf[]->MaterialContact() for each crack velocity
		//	field on that node
		nd[i]->MaterialContactOnNode(postUpdate,deltime);
	}
}

// adjust momenta at overlaping material velocity fields
void NodalPoint::CombineRigidMaterials(void)
{	// combine rigid materials across cracks
	int i;
	for(i=1;i<=nnodes;i++)
		nd[i]->CombineRigidParticles();
}


// Find Grid point velocities
void NodalPoint::GetGridVelocitiesForStrainUpdate(void)
{	int i;
    for(i=1;i<=nnodes;i++)
        nd[i]->CalcVelocityForStrainUpdate();
}

// Find Grid CM velocities (only for cracks when when contact.GetMoveOnlySurfaces() is FALSE)
void NodalPoint::GetGridCMVelocitiesTask8(void)
{	int i;
    for(i=1;i<=nnodes;i++)
        nd[i]->CalcCMVelocityTask8();
}


