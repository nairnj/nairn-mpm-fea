/********************************************************************************
	CrackVelocityFieldMulti.cpp
	nairn-mpm-fea

	Created by John Nairn on 21 August 2009.
	Copyright (c) 2009 John A. Nairn, All rights reserved.
********************************************************************************/
#if defined ( _MSC_VER) || defined (__APPLE__) 
#include "stdafx.h"
#endif
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Nodes/CrackVelocityFieldMulti.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Materials/ContactLaw.hpp"
#include "Nodes/MaterialContactNode.hpp"
#include "Exceptions/MPMWarnings.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "Elements/ElementBase.hpp"
#include <cmath>
#include <iostream>
using namespace std;

#pragma mark INITIALIZATION

// constructor
CrackVelocityFieldMulti::CrackVelocityFieldMulti(int num,short theLoc,int cnum) : CrackVelocityField(num,theLoc,cnum)
{	numberMaterials=0;
	numberRigidPoints=0;
}

// Destructor
CrackVelocityFieldMulti::~CrackVelocityFieldMulti()
{	// delete any allocated material velocity fields
	// go backwards so unmirror fields before source gets deleted
	for(int i=maxMaterialFields-1;i>=0;i--)
    {   if(mvf[i]!=NULL)
		{   // Always delete field [0], delete [i>0] if not ignoring cracks
			if(MVFInMemory(i))
				delete mvf[i];
			else
				mvf[i] = NULL;
		}
	}
}

// zero all active material velocity fields
// but mirrored fields just set to NULL
void CrackVelocityFieldMulti::ZeroMatFields(void)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(mvf[i]!=NULL)
		{	// If mirrored field that ignores cracks, just set to NULL
			if(!MVFInMemory(i))
				mvf[i] = NULL;
			else if(mvf[i]->numberPoints>0)
				mvf[i]->Zero();
		}
	}
	numberMaterials=0;
	numberRigidPoints=0;
}

// match material velocity fields on ghost node to those on a real node
// WARNING: this methods assumes no mirrored fields, which should be true because they are
//    not created until mass and momentem task and they are unmirrored at start of
//    initialization phase
// throws std::bad_alloc
void CrackVelocityFieldMulti::MatchMatVelocityFields(MatVelocityField **rmvf)
{	for(int i=0;i<maxMaterialFields;i++)
	{	if(rmvf[i]==NULL) continue;
        
		if(mvf[i]==NULL)
		{	mvf[i] = CreateMatVelocityField(rmvf[i]->GetFlags());
		}
		
		mvf[i]->Zero();
	}
	numberMaterials=0;
	numberRigidPoints=0;
}

#pragma mark TASK 1 AND 6 METHODS

// Called in intitation to preallocate material velocituy fields
// throws std::bad_alloc
void CrackVelocityFieldMulti::AddMatVelocityField(int matfld)
{	if(mvf[matfld]==NULL)
	{   mvf[matfld] = CreateMatVelocityField(MaterialBase::GetMVFFlags(matfld));
	}
}

// to avoid too much critical code, see if need to add mat velocity field first
bool CrackVelocityFieldMulti::NeedsMatVelocityField(int matfld) const { return mvf[matfld]==NULL; }

// add "mass" for rigid particle (task 1) - it counts particles too
void CrackVelocityFieldMulti::AddMassTask1(int matfld,double mnode,int numPts)
{	mvf[matfld]->mass += mnode;
	numberRigidPoints += numPts;
}

// Add to volume gradient. Do not call unless mpmgrid.extrapolateVolumeGradient>=0
// This gradient is only used in multimaterial contact calculations
// Corrections for symmetry planes and axisymmetry are done in GetVolumeGradient(). It is
//		faster there because that is called less while this is called for every node-particle pair
void CrackVelocityFieldMulti::AddVolumeGradient(int matfld,MPMBase *mptr,double dNdx,double dNdy,double dNdz)
{	double Vp = mptr->GetVolume(DEFORMED_AREA_FOR_GRADIENT);
	Vector dN = MakeVector(dNdx,dNdy,dNdz);
	mvf[matfld]->AddContactVector(mpmgrid.volumeGradientIndex,&dN,Vp);
}

// call this when copy from ghost to real node to sume gradients
// Do not call unless mpmgrid.extrapolateVolumeGradient>=0
void CrackVelocityFieldMulti::AddVolumeGradient(int matfld,Vector *grad)
{	mvf[matfld]->AddContactVector(mpmgrid.volumeGradientIndex,grad);
}

// This is called for crack field [i]. Each material field in crack source [0] that ignores cracks should be copied
//    to the corresponding field in this crack velocity field
// For NairnMPM, a new field is created if needed (instead of just mirroring it) and then filled with copy
void CrackVelocityFieldMulti::MirrorFieldsThatIgnoreCracks(CrackVelocityFieldMulti *cvfSource)
{
	// get pointer to material velocity field in the target crack velocity field
	MatVelocityField **mvfSource = cvfSource->GetMaterialVelocityFields();
	
	// loop over velocity fields in the source (which is in field [0])
	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvfSource[i]))
		{	if(mvfSource[i]->IgnoresCracks())
			{	// copy just the pointer
				mvf[i] = mvfSource[i];
				int sourcePoints = mvfSource[i]->numberPoints;
			
				// add number of rigid points to this crack velocity field
				if(mvfSource[i]->IsRigidField())
					numberRigidPoints += sourcePoints;
				
				// add total points in this crack velocity field
				numberPoints += sourcePoints;
			}
		}
	}
}

// Copy mass and momentum from ghost node to real node
void CrackVelocityFieldMulti::CopyMassAndMomentum(NodalPoint *real)
{	for(int matfld=0;matfld<maxMaterialFields;matfld++)
	{	if(mvf[matfld]!=NULL)
			mvf[matfld]->CopyMassAndMomentum(real,fieldNum,matfld);
	}
}

// zero momentum and displacement at a node for new calculations
// but can do the calculation for rigid particles here
void CrackVelocityFieldMulti::RezeroNodeTask6(double deltaTime)
{	int i;
    for(i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveField(mvf[i]))
        {	if(!mvf[i]->IsRigidField())
			{	mvf[i]->RezeroNodeTask6();
            }
            else if(MVFInMemory(i))
            {   // for rigid particles, keep initial pk
                // can project displacement using current velocity because
                // particle mass is its volume
                // dnew = Sum (Vp*fpi*(d + v dt)) = dold + Sum (Vp*fpi*v*dt) = dold + pk*dt
				if(mpmgrid.displacementIndex>=0)
					mvf[i]->AddContactVector(mpmgrid.displacementIndex,&mvf[i]->pk,deltaTime);
				if(mpmgrid.positionIndex>=0)
					mvf[i]->AddContactVector(mpmgrid.positionIndex,&mvf[i]->pk,deltaTime);
           }
        }
    }
}

// Copy mass and momentum from ghost node to real node
void CrackVelocityFieldMulti::CopyMassAndMomentumLast(NodalPoint *real)
{	for(int matfld=0;matfld<maxMaterialFields;matfld++)
    {	if(mvf[matfld]!=NULL)
            mvf[matfld]->CopyMassAndMomentumLast(real,fieldNum,matfld);
    }
}

#pragma mark TASK 3 METHODS

// Add to force spread out over the materials so each has same extra accerations = f/M
// Only called by AddTractionForce() when cracks are an imperfect interface
// Only adds force to fields that see cracks
void CrackVelocityFieldMulti::AddFtotSpreadTask3(Vector *f)
{	int i;
	
	// special case for only one material
	if(numberMaterials==1)
	{	for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			{	mvf[i]->AddFtot(f);
				break;
			}
		}
	}
	
	// more than one material, add to nonrigid materials only
	// when some material ignore cracks, do not add to those materials
	else
	{	double totMass=GetTotalMass(true);
		for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveNonrigidSeesCracksField(mvf[i],true))
				mvf[i]->AddFtotScaled(f,mvf[i]->mass/totMass);
		}
	}
}

// Add gravity and body force at a node in g mm/sec^2 (non rigid only)
// For materials that ignore cracks, only add to source MVF in field [0]
void CrackVelocityFieldMulti::AddGravityAndBodyForceTask3(Vector *gridBodyForce)
{	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidSourceField(mvf[i],fieldNum))
			mvf[i]->AddGravityAndBodyForceTask3(gridBodyForce);
	}
}

#ifdef RESTART_OPTION
// check for velocity and acceleration are too high
bool CrackVelocityFieldMulti::IsTravelTooMuch(double dt,double maxDist) const
{   for(int i=0;i<maxMaterialFields;i++)
    {   if(MatVelocityField::ActiveNonrigidSourceField(mvf[i],fieldNum))
        {   if(mvf[i]->IsTravelTooMuch(dt,maxDist))
                return true;
        }
    }
    return false;
}
#endif

// Copy grid forces ghost node to the real node (nonrigid only)
void CrackVelocityFieldMulti::CopyGridForces(NodalPoint *real)
{	for(int matfld=0;matfld<maxMaterialFields;matfld++)
	{	if(mvf[matfld]!=NULL)
            mvf[matfld]->CopyGridForces(real,fieldNum,matfld);
	}
}

// Restore momenta just prior to momentum update and to setting forces
// for velocity BCs
void CrackVelocityFieldMulti::RestoreMomenta(void)
{	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	mvf[i]->RestoreMomenta();
		}
	}
}

#pragma mark TASK 4 METHODS

// update momenta for this MPM step
//  pk(i+1) = pk(i) + ftot * dt
void CrackVelocityFieldMulti::UpdateMomentum(double timestep)
{	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidSourceField(mvf[i],fieldNum))
			mvf[i]->UpdateMomentum(timestep);
    }
}

#pragma mark XPIC METHODS

// Support XPIC calculations
void CrackVelocityFieldMulti::XPICSupport(int xpicCalculation,int xpicOption,NodalPoint *real,double timestep,int m,int k,double vsign)
{	if(xpicCalculation==COPY_VSTARNEXT) xpicOption = fieldNum;	// vfld
	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidSourceField(mvf[i],fieldNum))
		{	if(xpicCalculation==COPY_VSTARNEXT) m = i;		// matfld
			mvf[i]->XPICSupport(xpicCalculation,xpicOption,real,timestep,m,k,vsign);
		}
	}
}

#pragma mark MATERIAL CONTACT

/* Called in multimaterial mode to check contact at nodes with multiple materials

	Input parameters:
		material mass, displacement and/or position, volume gradient (if needed)

	Output changes are only allowed on this node (to be thread safe for parallel)
		changes on mvf[]: pk, ftot (if postUpdate is TRUE)
 
	ndptr is parent node to this crack velocity field
 
	On first call in time step, first and last are pointers to Cracknode * because those
		objects are created for later interface calculations
	postUpdate is TRUE when called between momentum update and particle update and otherwise is false
	throws std::bad_alloc
*/
void CrackVelocityFieldMulti::MaterialContactOnCVF(MaterialContactNode *mcn,double deltime,int callType)
{
	// exit if no contact
	if(numberMaterials<=1) return;
	
	// list of active non-rigid materials
	int *activeMat = new int[numberMaterials];
	int numMats = 0;
	
	// get center of mass results and look out for rigid materials
	// uses nodal mass and mommentum
	int i,rigidMat=-1;
	bool multiRigid = false;
	Vector Pc;
	ZeroVector(&Pc);
	double Mc=0.;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	if(!mvf[i]->IsRigidField())
			{	// real material
				AddVector(&Pc,&mvf[i]->pk);
				Mc += mvf[i]->mass;
				activeMat[numMats++] = i;
			}
			else if(rigidMat>=0)
			{	// rigid material, but not allowed if already had another rigid material
				multiRigid = true;
			}
			else
			{	// first rigid material at this node
				rigidMat=i;
			}
		}
	}
	
	// if rigid materials present, then special case for contact laws and then done
	if(rigidMat>=0)
	{	RigidMaterialContactOnCVF(rigidMat,multiRigid,mcn,deltime,callType);
	}
	else
	{
#ifdef THREE_MAT_CONTACT
		// Must specify lumping methods in multimaterial command (only in XML for now)
		if(numMats>2 && mpmgrid.lumpingMethod==EXPLICIT_PAIRS)
		{	int handled = MaterialContactThree(mcn,deltime,callType,activeMat,numMats);
			if(DisplayThreeMatWarnings(handled))
				MaterialContactOnCVFLumped(mcn,deltime,callType,activeMat,numMats,Pc,Mc);
		}
		else
			MaterialContactOnCVFLumped(mcn,deltime,callType,activeMat,numMats,Pc,Mc);
#else
		MaterialContactOnCVFLumped(mcn,deltime,callType,activeMat,numMats,Pc,Mc);
#endif // end THREE_MAT_CONTACT
	}
	
	delete [] activeMat;
	return;
}

// Contact for each material lumped with the other materials
void CrackVelocityFieldMulti::MaterialContactOnCVFLumped(MaterialContactNode *mcn,double deltime,
									int callType,int *activeMat,int numMats,Vector Pc,double Mc)
{
	// get the node
	NodalPoint *ndptr = mcn->GetTheNode();
	
	double qrate = 0.;
	int i;
	
	// loop materials or pairs of materials
	bool doingPairs = (numMats==2 && mpmgrid.materialNormalMethod!=EACH_MATERIALS_MASS_GRADIENT);	// has only one pair;
	int miMax = doingPairs ? numMats-1 : numMats ;
	
	// Get center of mass displacememnt used in this mode
	Vector dispc = GetCMDisplacement(ndptr,false,mpmgrid.contactByDisplacements);
	ScaleVector(&dispc,1./Mc);
	Vector dispcForInterface = dispc;
	if(!mpmgrid.contactByDisplacements && mpmgrid.hasImperfectInterface)
	{	// get displacement-based result for interfaces
		dispcForInterface = GetCMDisplacement(ndptr,false,true);
		ScaleVector(&dispcForInterface,1./Mc);
	}
	
	MatVelocityField *mvfi,*mvfj;
	for(int mi=0;mi<miMax;mi++)
	{	// known active field - get volume from mvf->volume
		i = activeMat[mi];
		mvfi = mvf[i];
		double massi = mvfi->mass;
		double voli = mvfi->GetContactVolume();
		double volj = GetContactVolumeNonrigid(false) - voli;

		// Contact law from other material with most volume (uses mvf->volume)
		int j = -1;
		
		// If needed, find volume gradient of all other non-rigid materials
		Vector gradj;
		if(mpmgrid.volumeGradientIndex>=0) ZeroVector(&gradj);
		double maxOtherMaterialVolume = 0.;
		for(int kj=0;kj<numMats;kj++)
		{	if(kj==mi) continue;
			int jj = activeMat[kj];
			double matVolume = mvf[jj]->GetContactVolume();
			if(matVolume>maxOtherMaterialVolume)
			{	maxOtherMaterialVolume = matVolume;
				j = jj;
			}
			// Finding  - Sum grad V_j
			if(mpmgrid.volumeGradientIndex>=0)
			{	Vector normj;
				GetVolumeGradient(jj,ndptr,&normj,-1.);
				AddVector(&gradj,&normj);
			}
		}
		
		// problem if paired material not found, but it will always be found
		if(j<0) continue;
		mvfj = mvf[j];
		ContactLaw *theContactLaw = mpmgrid.GetMaterialContactLaw(i,j);
		
		// compared to lumped mass (or other material)
		double mred = (Mc-massi)/Mc;		// final mred = massi*(Mc-massi)/Mc, but found later when needed
		double massRatio = massi/Mc;
		
		// some variables
		Vector norm,delta;					// normal and COD (later if needed)
		double dotn,deln = 0.;				// normal force and cod components
		Vector delMats,*delMatsPtr=NULL;	// for distances to materials in get normal call
		double contactArea = -1.,mredDelWf = -1.,contactGridN = 1.;
		
		// Determine if in contact
		bool comContact = false;
		
		// find mi(vc-vi) = (mi/mc)pc-pi = (mi pc - mc pi)/mc or momentum change to match ctr of mass momentum
		// second version might be less prone to round-off error if low or high massRatio
		Vector delPi;
		CopyScaleVector(&delPi,&mvfi->pk,-1.);
		AddScaledVector(&delPi,&Pc,massRatio);
		//CopyScaleVector(&delPi,&Pc,massi);
		//AddScaledVector(&delPi,&mvfi->pk,-Mc);
		//ScaleVector(&delPi,1./Mc);
		AdjustForSymmetry(ndptr, &delPi, false);
		
		// If needed, determine if the surfaces are in contact
		if(theContactLaw->IgnoreContact() || comContact)
		{	// always in contact and use center of mass to stick
			comContact = true;
		}
		
		else
		{	// Ignore very small mass nodes - may not be needed
			if(massRatio<MASS_MIN || massRatio>MASS_MAX)
			{	continue;
			}

			// Get normal vector by various options
			bool hasDeln = false;
			if(!NonRigidCustomNormal(ndptr,i,j,norm))
			{	norm = GetNormalVector(mcn,i,-1,voli,&gradj,volj,1.,hasDeln,&delMats);
				if(hasDeln)
				{	deln = delMats.x;
					delMatsPtr = &delMats;
				}
				
				// nan here means volume gradients zero as gradient normal to a symmetry plane
				if(norm.x!=norm.x || norm.y!=norm.y || norm.z!=norm.z) continue;
			}
			
			// Get dotn = sticking force (actual (vc-vi) = delPi/mi and (vb-va) = delPi/mred (when mred is done)
			dotn=DotVectors(&delPi,&norm);
			
			// Get deln (if not done by LR). Also get delta and contactArea if needed
			// Note that all interface need contactArea; they use delta later to get tangDel
			if(!hasDeln || theContactLaw->ContactLawNeedsContactArea())
			{	// Get displacement for material a
				Vector dispa;
				CopyScaleVector(&dispa,mvfi->GetContactDispPtr(mpmgrid.contactByDisplacements),1./massi);
				AdjustForSymmetry(ndptr,&dispa,false);
				
				// Get COD based on other material lumped
				// Let dispc = xiCOM, dispa = xia = (displacement or position)/massi, and mred* = (Mc-Mi)/Mc
				// ... then from my 2013 paper, COD = xib - xia = (Mc/(Mc-mi))(xiCOM - xia) = (dispc-xia)/mred*
				delta = dispc;
				SubVector(&delta,&dispa);
				ScaleVector(&delta,1./mred);
				
				// get deln = delta.n and correct it if simulation is extrapolating positions
				if(!hasDeln)
				{	// get dispb = delta + xia = xib-xia+xia = xib (already adjusted for symmetry)
					Vector dispb = delta;
					AddVector(&dispb,&dispa);
					deln = contact.MaterialSeparation(DotVectors(&dispb,&norm),DotVectors(&dispa,&norm),&norm,ndptr,
													  mpmgrid.contactByDisplacements,mpmgrid.positionCutoff);
				}
				
				// for interfaces or contact laws that need, get area here
				if(theContactLaw->ContactLawNeedsContactArea())
				{	contactArea = GetContactArea(ndptr,voli,volj,&norm,&contactGridN,delMatsPtr);
					
					if(!mpmgrid.contactByDisplacements && theContactLaw->IsImperfectInterface())
					{	// recalculate delta using displacements if needed for an interface
						CopyScaleVector(&dispa,mvfi->GetContactDispPtr(true),1./massi);
						AdjustForSymmetry(ndptr,&dispa,false);
						delta = dispcForInterface;
						SubVector(&delta,&dispa);
						ScaleVector(&delta,1./mred);
					}
				}
			}
			
			// Final reduced mass = (mi*(Mc-mi))/Mc while until here was (Mc-mi)/Mc
			mred *= massi;
        }
			
		// Here inContact is usually true, but may be false for imperfet interface or for
		// any other contact law that want to continue (e.g., friction with adhesion)
		if(comContact)
		{	// Use delPi found above and no further change needed, no other data needed either
		}
		else if(theContactLaw->IsFrictionalContact())
		{	// Handle frictional law
			
			// Find dF = (ma Fb - mb Fa)/Mc = (ma Fc - Mc Fa)/Mc = (ma Fc/Mc) - Fa (when needed)
			// note that dacc = dF/mred
			bool getHeating = false;
			Vector delFi;
			Vector *delFiPtr = NULL;
			if(callType==UPDATE_MOMENTUM_CALL)
			{	delFi = GetCMatFtot();
				ScaleVector(&delFi, massi/Mc);
				AddScaledVector(&delFi, mvfi->GetFtotPtr(), -1.);
				delFiPtr = &delFi;
				// second order heating needs acceleration too
				if(ConductionTask::matContactHeating) getHeating = true;
			}
			
			// get change in momentum, but false return means friction law want no change (i.e, continue as if not in contact)
			if(!theContactLaw->GetFrictionalDeltaMomentum(&delPi,&norm,dotn,deln,&mredDelWf,mred,getHeating,
														  contactArea,deltime,delFiPtr,ndptr))

			{	continue;
			}
		}
		else
		{	// Handle imperfect interface
			
			// tangential COD
			Vector tangDel;
			double delt = GetTangentCOD(&norm,&delta,&tangDel);
			
			// Find delFi = (ma Fb - mb Fa)/Mc = (ma Fc - Mc Fa)/Mc = (ma Fc/Mc) - Fa (when needed)
			Vector fImp;
			if(callType==UPDATE_MOMENTUM_CALL)
			{	fImp = GetCMatFtot();
				ScaleVector(&fImp, massi/Mc);
				AddScaledVector(&fImp, mvfi->GetFtotPtr(), -1.);
			}
			else
				ZeroVector(&fImp);
			
			// Get interface force, energy, and altered delPi
			double rawEnergy;
			theContactLaw->GetInterfaceForces(&norm,&fImp,&rawEnergy,contactArea,&delPi,dotn,mred,
											  &tangDel,deln,delt,contactGridN);
			if(callType==UPDATE_MOMENTUM_CALL)
			{	// add force (if any) to momentum change
				AddScaledVector(&delPi, &fImp, timestep);
				
				// Add interface energy. (Legacy units g-mm^2/sec^2 or multiply by 1e-9 to get J - kg-m^2/sec^2)
#pragma omp atomic
				NodalPoint::interfaceEnergy += rawEnergy;
			}
		}
		
#ifdef CHECK_NAN
		if(delPi.x!=delPi.x || delPi.y!=delPi.y || delPi.z!=delPi.z)
		{
#pragma omp critical (output)
			{	cout << "\n# CrackVelocityFieldMulti::MaterialContactOnCVF: bad delPi";
				PrintVector(" = ",&delPi);
				PrintVector(" norm = ",&norm);
				cout << endl;
				ndptr->Describe(true);
			}
		}
#endif
		
		// Apply momentum change to the velocity fields

		// change momenta
		mvfi->ChangeMatMomentum(&delPi,callType,deltime);
		
		// special case two materials for efficiency (and if both will find normal the same way)
		if(doingPairs)
		{	mvfj->ChangeMatMomentum(ScaleVector(&delPi,-1.),callType,deltime);
			
			// only true when in momentum update and conduction is on
			if(mredDelWf>0.)
			{	qrate += mredDelWf/mred;
			}
		}
		else
		{	// separate lumped materials
		
			// friction when more than two materials (only true when in momentum update and conduction is on)
			if(mredDelWf>0.)
			{   qrate += mredDelWf/massi;
			}
		}

	}
	
	// total friction, scale and absolute value
	if(callType==UPDATE_MOMENTUM_CALL && ConductionTask::matContactHeating && qrate!=0.)
	{	// add to flux
		conduction->AddFluxCondition(ndptr,fabs(qrate/deltime),true);
		
#pragma omp atomic
		NodalPoint::frictionWork += qrate;
	}
}

// Called in multimaterial mode to check contact at nodes with multiple materials and here
// contact with one or more rigid materials (no rigid materials handled in MaterialContactOnCVF()
// throws std::bad_alloc
void CrackVelocityFieldMulti::RigidMaterialContactOnCVF(int rigidFld,bool multiRigid,
										MaterialContactNode *mcn,double deltime,int callType)
{
	// multiRigid means more than one rigid material
	int i;
	double rigidVolume = mvf[rigidFld]->GetContactVolume();
	if(multiRigid)
	{   // this aborts the calculation
		//throw CommonException("Two different rigid materials in contact on the same node",
		//						"CrackVelocityFieldMulti::RigidMaterialContactOnCVF");
		
		// scheme 1, find material with the most contact volume
		for(i=rigidFld+1;i<maxMaterialFields;i++)
		{   if(MatVelocityField::ActiveField(mvf[i]))
			{	if(mvf[i]->IsRigidField())
				{   double testVolume = mvf[i]->GetContactVolume();
					if(testVolume > rigidVolume)
					{   rigidVolume = testVolume;
						rigidFld = i;
					}
				}
			}
		}
		
		// alternate scheme might combine the rigid particles in some way, but a problem
		// with two rigid particles approaching each other could end up with zero average
		// normal vector and therefore poor results.
	}
	
	// get rigid material volume for contact and actual volume (rigidVolume is above)
	// these are the same except in axisymmtric where first is area and second is volume per radian
	double actualRigidVolume = mvf[rigidFld]->mass;
	
	// rigid material velocity with position-dependent velocities may have mixed velocities
	Vector rvel;
	CopyScaleVector(&rvel, &mvf[rigidFld]->pk, 1./actualRigidVolume);
	
	// loop over each material (skipping the rigid materials)
	NodalPoint *ndptr = mcn->GetTheNode();
	for(i=0;i<maxMaterialFields;i++)
	{	// skip if inactive or if a rigid material
		if(!MatVelocityField::ActiveField(mvf[i]) || mvf[i]->IsRigidField()) continue;
		
		// First determine contact law with rigid material
		ContactLaw *theContactLaw = mpmgrid.GetMaterialContactLaw(i,rigidFld);
		
		// get contact parameters
		Vector norm,delta;				// normal and COD (if needed)
		Vector dispRigid;				// displacement of rigid material (might be needed for Yang correction)
		double dotn,deln=0.;			// normal force and cod
		bool hasDeln;					// might be needed for Yang correction
		Vector delMats,*delMatsPtr=NULL;		// for distances to the node in get normal calls
		double contactArea = -1.,mredDelWf = -1.,contactGridN = 1.;
 	
		// some variables
		double massi=mvf[i]->mass;
		double mred=massi;
		double voli=mvf[i]->GetContactVolume();
		
		// momentum change
		Vector delPi;
		// find mi(vr-vi) = - pi + mi*vr or momentum change to match rigid particle velocity
        CopyScaleVector(&delPi,&mvf[i]->pk,-1.);
        AddScaledVector(&delPi,&rvel,massi);
		AdjustForSymmetry(ndptr,&delPi,false);
		
		// if needed, determine if the surfaces are in contact
		bool comContact = false;
		if(theContactLaw->IgnoreContact())
		{   // will use delPi to match rigid particle velocity
			comContact = true;
		}
		
		else
		{   // Ignore very small interactions
			double volRatio=voli/(rigidVolume+GetContactVolumeNonrigid(false));
			if(volRatio<MASS_MIN || volRatio>MASS_MAX)
			{	continue;
			}
			
			// Get normal vector, dotn, and deln
			// Also find displacment vector
			
			// Get normal vector by various options
			hasDeln = false;
			if(!RigidCustomNormal(ndptr,i,rigidFld,norm))
			{	Vector gradj;
				if(mpmgrid.volumeGradientIndex>=0)
					GetVolumeGradient(rigidFld,ndptr,&gradj,-1.);
				norm = GetNormalVector(mcn,i,rigidFld,voli,&gradj,rigidVolume,mpmgrid.rigidGradientBias,hasDeln,&delMats);
				if(hasDeln)
				{	deln = delMats.x;
					delMatsPtr = &delMats;
				}
				
				// nan here means volume gradients zero as gradient normal to a symmetry plane
				if(norm.x!=norm.x || norm.y!=norm.y || norm.z!=norm.z) continue;
			}
			
			// get approach direction momentum from delPi.n (actual (vc-vi) = delPi/mi)
			dotn=DotVectors(&delPi,&norm);
			
			// Get deln (if not done by LR). Also get delta and contactArea if needed
			// Note that all interface need contactArea; they use delta later to get tangDel
			if(!hasDeln || theContactLaw->ContactLawNeedsContactArea())
			{	// Displacement calculations: rigid material displacement was scaled by volume, while non-rigid was weighted by mass
				CopyScaleVector(&dispRigid,mvf[rigidFld]->GetContactDispPtr(mpmgrid.contactByDisplacements),1./actualRigidVolume);
				AdjustForSymmetry(ndptr,&dispRigid,false);
				Vector dispi;
				CopyScaleVector(&dispi,mvf[i]->GetContactDispPtr(mpmgrid.contactByDisplacements),1./massi);
				AdjustForSymmetry(ndptr,&dispi,false);
				
				// get deln = delta.n corrected if simulation is extrapolating positions or displacements
				if(!hasDeln)
				{	deln = contact.MaterialSeparation(DotVectors(&dispRigid,&norm),DotVectors(&dispi,&norm),&norm,ndptr,
													  mpmgrid.contactByDisplacements,mpmgrid.positionCutoff);
				}
				
				if(theContactLaw->ContactLawNeedsContactArea())
				{	// interface and some laws will need contact area
					contactArea = GetContactArea(ndptr,voli,rigidVolume,&norm,&contactGridN,delMatsPtr);
					
					// get COD if might be nneded by an interface
					if(!mpmgrid.contactByDisplacements && mpmgrid.hasImperfectInterface)
					{	CopyScaleVector(&dispRigid,mvf[rigidFld]->GetContactDispPtr(true),1./actualRigidVolume);
						AdjustForSymmetry(ndptr,&dispRigid,false);
						CopyScaleVector(&dispi,mvf[i]->GetContactDispPtr(true),1./massi);
						AdjustForSymmetry(ndptr,&dispi,false);
					}
					delta = dispRigid;
					SubVector(&delta,&dispi);
				}
			}
		}
		
		// Here inContact is usually true, but may be false for imperfect interface or for
		// any other contact law that wants to continue (e.g., friction with adhesion)
		if(comContact)
		{	// Use delPi found above and no further change needed, no other data needed either
		}
		else if(theContactLaw->IsFrictionalContact())
		{	// Handle frictional law
			
			bool getHeating = false;
			Vector delFi;
			Vector *delFiPtr = NULL;
			if(callType==UPDATE_MOMENTUM_CALL)
			{	// Find dF = (ma Fb - mb Fa)/Mc = -Fa for rigid material (infinite mb)
				// May not be correct for rigid block, but not sure what else to do
				// note that dacc = dF/mred = -Fa/ma
				delFi = mvf[i]->GetFtot();
				ScaleVector(&delFi, -1.);
				delFiPtr = &delFi;
				// second order heating needs acceleration too
				if(ConductionTask::matContactHeating) getHeating = true;
			}
			
			// get change in momentum, but false return means friction law wants no change (i.e, continue as if not in contact)
			if(!theContactLaw->GetFrictionalDeltaMomentum(&delPi,&norm,dotn,deln,&mredDelWf,mred,getHeating,
														  contactArea,deltime,delFiPtr,ndptr))
			{	continue;
			}

//#define YANG_CORRECTION
//#define USE_DELN_INSTEAD
#ifdef YANG_CORRECTION
			// To fully implement have commands to set this parameters with enum for decayFunction
			// options decayFunction=0 turns it off and > last function reverts to truncation at contrainStart
			double constrainStart = 0;
			int decayFunction = 6;
			double kDecay = 8.;					// sigmoidal seems worse than linear so maybe not needed
			
			if(decayFunction>0)
			{	// The following section is inspired by Yang theses (U of W with Peter Arduino) which scaled
				//    force (i.e. delPi) depending on distance from node to edge of the rigid material
				//    Yang defined a rigid body. Here rigid body is collection of rigid particles
				// It scales delPi by a factor from 0 to 1 as function of distance of the rigid particle
				//    to the nodal point. It might be an approach to soften the contact phenomenon,
				//    or made it should be softened based on distance between the materials?
				// Here modified to various functions and allow full contraint cutoff to differ from zero
				// Wang thesis recommends contrainStart=0 and decayFunction=2 (or (1-r/h)^2)
				// Shock examples are improved with contrainStart=-1 and decayFunction=1 (or 1-r/h)
#ifdef USE_DELN_INSTEAD
				double db = deln/contactGridN;
				
				// decay function = 1 (db<constrainStart, which must be negative), f(db) (db<0), 0 (db>=0)
				if(db>=0.)
				{	// not in contact so no contact force
					ZeroVector(&delPi);
				}
#else
				double db;
				if(delMatsPtr!=NULL)
				{	// using logistic regression
					db = delMats.z/contactGridN;
				}
				else if(mpmgrid.contactByDisplacements)
				{	// no correction here because don't have actual separation
					db = constrainStart-1.;
				}
				else
				{	// find extrapolated distance from rigid material to the node using contactPosition setting
					// Rigid material is material b snf (a,b) contact pair
					double pb = (dispRigid.x-ndptr->x)*norm.x + (dispRigid.y-ndptr->y)*norm.y + (dispRigid.z-ndptr->z)*norm.z;
					pb /= contactGridN;
					double r = mpmgrid.positionCutoff;
					if(r>0.)
					{	// contact position in units of cells
						db = pb - 0.5*r;
					}
					else
					{	// fit to non-linear power law
						r = -r;
						db = 2.*pow(pb/1.25,r) - 1.;
					}
				}
				
				// decay function = 1 (db<constrainStart), f(db) (db<1), 0 (db>1)
				if(db>1.)
				{	// more than grid cell away so no contact force
					ZeroVector(&delPi);
				}
#endif
				else if(db>constrainStart)
				{	// rescale for step function (r/h = 0 when db=constrainStart and 1 when db=0)
#ifdef USE_DELN_INSTEAD
					double roverh = 1. - db/constrainStart;
#else
					double roverh = (db-constrainStart)/(1.-constrainStart);
#endif
					// various decay functions
					double decay;
					switch(decayFunction)
					{	case 1://LINEAR_DECAY
							decay = 1.-roverh;
							break;
						case 2://SQUARE_DECAY
							// advocated in Yang's thesis
							decay = (1.-roverh)*(1.-roverh);
							break;
						case 3://SQUARE2_DECAY
							decay = 1.-roverh*roverh;
							break;
						case 4://CUBICSMOOTHSTEP_DECAY
							decay = 1.-roverh*roverh*(3.-2.*roverh);
							break;
						case 5://Sigmoidal
						{	double x = 2*roverh-1.;			// rescale to -1 to 1
							double sx = 1./(1+exp(kDecay*x));
							double sm1 = 1./(1+exp(-kDecay));
							double s1 = 1./(1+exp(kDecay));
							decay = (sx-s1)/(sm1-s1);
							break;
						}
						case 6://Conserve constant velocity gradient for one side (set contrainStart to 0)
							decay = pow((1.-roverh),2.5);
							break;
						default: //TRUNCATE_DECAY
							decay = 0.;
							break;
					}
//#pragma omp critical (debug)
					//{	cout << "deln=" << deln << ", db=" << db << ", r/h=" << roverh << ", decay=" << decay << endl;
					//}

					// scale normal contact force only
					double dpnorm = DotVectors(&delPi,&norm);
					AddScaledVector(&delPi,&norm,(decay-1.)*dpnorm);
				}
				// else db<constrainStart and no change or fully constrained
			}
#endif
		}
		else
		{	// Handle imperfect interface
			Vector tangDel;
			
			// tangential COD
			double delt = GetTangentCOD(&norm,&delta,&tangDel);
			
			Vector fImp;
			ZeroVector(&fImp);
			if(callType==UPDATE_MOMENTUM_CALL)
			{	// Find delFi = (ma Fb - mb Fa)/Mc = (ma Fc - Mc Fa)/Mc = (ma Fc/Mc) - Fa = -Fa (when needed)
				// not sure best for rigid black, but don't know what else to do
				AddScaledVector(&fImp, mvf[i]->GetFtotPtr(), -1.);
			}
			
			// Get interface force, energy, and posible change delPi
			double rawEnergy;
			theContactLaw->GetInterfaceForces(&norm,&fImp,&rawEnergy,
													contactArea,&delPi,dotn,mred,&tangDel,deln,delt,contactGridN);
			if(callType==UPDATE_MOMENTUM_CALL)
			{	// add force (if any) to momentum change
				AddScaledVector(&delPi, &fImp, timestep);
				
				// Add interface energy. (Legacy units g-mm^2/sec^2 or multiply by 1e-9 to get J - kg-m^2/sec^2)
#pragma omp atomic
				NodalPoint::interfaceEnergy += rawEnergy;
			}
		}
		
#ifdef CHECK_NAN
		if(delPi.x!=delPi.x || delPi.y!=delPi.y || delPi.z!=delPi.z)
		{
#pragma omp critical (output)
			{	cout << "\n# CrackVelocityFieldMulti::MaterialContactOnCVF: bad delPi";
				PrintVector(" = ",&delPi);
				PrintVector(" norm = ",&norm);
				cout << endl;
				ndptr->Describe(true);
			}
		}
#endif
		
		// Change momentum on field i
		mvf[i]->ChangeMatMomentum(&delPi,callType,deltime);
		
		// store contact force in rigid particle ftot
		// There are three times this is called in an MPM Step, but to get correct force,
		//    this calculation must only be done during momentum update
		if(callType==UPDATE_MOMENTUM_CALL)
		{	mvf[rigidFld]->AddContactForce(&delPi);
			
			// if conduction on (and here means in momentum update) add frictional heating
			if(mredDelWf>0. && ConductionTask::matContactHeating)
			{	// add to flux
				double qrate = mredDelWf/mred;
				conduction->AddFluxCondition(ndptr,fabs(qrate/deltime),true);
				
#pragma omp atomic
				NodalPoint::frictionWork += qrate;
			}
		}
	}
}

// Get normal vector from material i to j by various options
// Note 1: j<0 means others lumped. j>0 is for rigid contact, or explicit 3 (not recommended anymore)
// Note 2: gradj only defined when volumeGradientIndex>=0
// Note 3: If LR method, hasDeln will change to true and set delMats to
//			(dBeta-dAlpha,dAlpha,dBeta) for separation and signed distances of each material to node
Vector CrackVelocityFieldMulti::GetNormalVector(MaterialContactNode *mcn,int i,int j,double voli,Vector *gradj,double otherVolume,
												double jBias,bool &hasDeln,Vector *delMats)
{
	NodalPoint *ndptr = mcn->GetTheNode();
	Vector norm;
	
	// Pick and option
	switch(mpmgrid.materialNormalMethod)
	{	case MAXIMUM_VOLUME_GRADIENT:
		{	// Use mat with largest magnitude volume gradient
			Vector normi;
			GetVolumeGradient(i,ndptr,&normi,1.);
			
			// compare magnitude of volume gradients (looking at squares and jBias is squared)
			double magi = sqrt(DotVectors(&normi,&normi));
			double magj = sqrt(DotVectors(gradj,gradj));
			if(magi >= jBias*magj)
				CopyScaleVector(&norm,&normi,1./magi);			// use material i
			else
				CopyScaleVector(&norm,gradj,1./magj);			// lumped other material
			break;
		}
			
		case MAXIMUM_VOLUME:
			// Use mat with most volume
			if(voli >= otherVolume)
				GetVolumeGradient(i,ndptr,&norm,1.);
			else
				CopyVector(&norm,gradj);
			CopyScaleVector(&norm,&norm,1./sqrt(DotVectors(&norm,&norm)));
			break;
			
		case AVERAGE_MAT_VOLUME_GRADIENTS:
		{	// avg grad is (gradi + gradj)/(voli+otherVolume)
			GetVolumeGradient(i,ndptr,&norm,1.);
			AddVector(&norm,gradj);
			double magi = DotVectors(&norm,&norm);
			if(j>=0)
			{	if(mvf[j]->IsRigidField())
				{	// rigid grad is gradj/otherVolume
					double vRatio = (voli + otherVolume)/otherVolume;
					double magj = DotVectors(gradj,gradj);
				
					// compare square of volume gradients (jBias has been squared)
					if(magi >= jBias*magj*vRatio*vRatio)
						ScaleVector(&norm,1./sqrt(magi));				// use the average
					else
						CopyScaleVector(&norm,gradj,1./sqrt(magj));		// use rigid material
				}
				else
				{	// just normalize
					ScaleVector(&norm,1./sqrt(magi));
				}
			}
			else
			{	// just normalize
				ScaleVector(&norm,1./sqrt(magi));
			}
			break;
		}
			
		case EACH_MATERIALS_MASS_GRADIENT:
			// Use each material's own gradient and handle separately (i.e. does not conserve momentum)
			GetVolumeGradient(i,ndptr,&norm,1.);
			ScaleVector(&norm,1./sqrt(DotVectors(&norm,&norm)));
			break;
			
		case SPECIFIED_NORMAL:
			// use specified normal for all contact
			norm = mpmgrid.contactNormal;
			AdjustForSymmetry(ndptr,&norm,true);
			break;

		case LINEAR_REGRESSION:
			if(jBias>=100.)
			{	double magj = DotVectors(gradj,gradj);
				if(DbleEqual(magj,0.))
					norm = LinearRegressionNormal(mcn,fieldNum,i,j,hasDeln,delMats);
				else
				{	CopyScaleVector(&norm,gradj,1./sqrt(magj));		// use rigid material
					FindSepFromNormalAndPointCloud(&norm,mcn,fieldNum,i,j,delMats);
					hasDeln = true;
				}
			}
			else
				norm = LinearRegressionNormal(mcn,fieldNum,i,j,hasDeln,delMats);
			break;
			
		case LOGISTIC_REGRESSION:
			if(jBias>=100.)
			{	double magj = DotVectors(gradj,gradj);
				if(DbleEqual(magj,0.))
					norm = LogisticRegressionNormal(mcn,fieldNum,i,j,hasDeln,delMats);
				else
				{	CopyScaleVector(&norm,gradj,1./sqrt(magj));		// use rigid material
					FindSepFromNormalAndPointCloud(&norm,mcn,fieldNum,i,j,delMats);
					hasDeln = true;
				}
			}
			else
				norm = LogisticRegressionNormal(mcn,fieldNum,i,j,hasDeln,delMats);
			break;
			
		default:
			break;
	}
	
	// final result
	return norm;
}

// This method allows the code to return a custom normal for contact between non-rigid materials.
// Typically a custom normal will be selected using a developer flag (in dFlag[]).
// If not custom, return false
bool CrackVelocityFieldMulti::NonRigidCustomNormal(NodalPoint *ndptr,int i,int j,Vector &norm)
{
	// no custom normals available
	return false;
}

// This method allows the code to return a custom normal for contact with rigid materials.
// Typically a custom normal will be selected using a developer flag (in dFlag[]).
// If not custom, return false
bool CrackVelocityFieldMulti::RigidCustomNormal(NodalPoint *ndptr,int i,int rigidFld,Vector &norm)
{
	// Various developer flag hacks to change the normal vector in specific problems
		
	return false;
}

// On input of all, return vector for interfacial opening displacement and load dispa and dispb
Vector CrackVelocityFieldMulti::GetDisplacementVector(NodalPoint *ndptr,Vector *dispaEx,double massa,Vector *dispbEx,double massb,
													  bool hasDeln,double deln,Vector *norm,Vector *dispa,Vector *dispb,bool &hasDelta)
{	// now have it
	hasDelta = true;
	
	// get dispa for this material from xib - xia
	CopyScaleVector(dispa,dispaEx,1./massa);
	CopyScaleVector(dispb,dispbEx,1./massb);
	Vector delta = *dispb;
	SubVector(&delta,dispa);
	
	if(hasDeln)
	{	// accept tangential component of delta, but change normal to deln
		double extrapDeln = DotVectors(&delta,norm);
		AddScaledVector(&delta,norm,deln-extrapDeln);
	}
	
	// adjust and return
	AdjustForSymmetry(ndptr,&delta,false);
	return delta;
}

// Get tangDel and delt (only used in interface laws)
double CrackVelocityFieldMulti::GetTangentCOD(Vector *norm,Vector *delta,Vector *tangDel)
{
	// tangential vector in tangDel (might be zero vector if delta || to n)
	double deln = DotVectors(delta,norm);
	CopyVector(tangDel,delta);
	AddScaledVector(tangDel,norm,-deln);				// delta - deln (n) = delt (t)
	
	// by normalizing to positive delt, hat t always points in positive direction
	double delt=sqrt(DotVectors(tangDel,tangDel));
	if(!DbleEqual(delt,0.)) ScaleVector(tangDel,1./delt);
	
	// return tangDel in that vector and tangDel.t in double
	return delt;
}

// Get interfacial contact area and also find tangDel and delt (only used in interface laws)
// Input is ndptr, voli, volb, norm
// Ouput is hperp, and return contactArea
// If tangDel!=NULL also find tangDel and delt = tangDel.t from delta and norm (otherwase delta not even used)
// If delMats!=NULL, it is (deln, dAlpha, dBeta) for separation and distance of each material to node)
double CrackVelocityFieldMulti::GetContactArea(NodalPoint *ndptr,double voli,double volb,Vector *norm,double *hperp,Vector *delMats)
{
	// perpendicular distance to correct contact area and contact by positions
	Vector dist = mpmgrid.GetPerpendicularDistance(norm,ndptr);
	
	// Get raw surface area, it is divided by hperp to get contact area
	// multiple by position for axisymmetric
	double surfaceArea;
	if(delMats==NULL)
	{	// Scale voltot=voli+volb to voltot*sqrt(2*vmin/voltot) = sqrt(2*vmin*vtot)
		// dist weightings to allow for unequal element sizes
		// This is approach in my paper and estimates distance from interface to the node
		if(ElementBase::GetShapeFunctionOrder()==1)
			surfaceArea = sqrt(2.*fmin(voli*dist.y,volb*dist.z)*(voli+volb))/dist.x;
		else
		{	// spline functions
			double Vcell = voli+volb;
			
			// get linear d function and convert for splines
			double ds = 1. - sqrt(2.*fmin(voli*dist.y,volb*dist.z)/Vcell);
			double d = fmax(1.193*ds,-4.3248 + ds*(18.5999 - ds*(23.43 - 10.6304*ds)));
			
			if(d<0.5)
				surfaceArea = 0.75 - d*d;
			else
			{	double d2 = fmax(0.,3. - 2*d);
				surfaceArea = 0.125*d2*d2;
			}
			
			// Now contact area is Vcell*S(d)/hperp
			surfaceArea *= Vcell/dist.x;
		}
	}
	else
	{	// absolute value of each material to the node
		double dAlpha = fabs(delMats->y);
		double dBeta = fabs(delMats->z);
		//double Vcell = mpmgrid.GetCellVolume(ndptr);
		double Vcell = voli+volb;
		
		// Find material closest to node in cell units (d>0)
		double d = fmin(dAlpha,dBeta)/dist.x;
		
		// get linear or quadratic weighting factor = S(d)
		if(ElementBase::GetShapeFunctionOrder()==1)
			surfaceArea = fmax(0.,1.-d);
		else
		{	if(d<0.5)
				surfaceArea = 0.75 - d*d;
			else
			{	double d2 = fmax(0.,3. - 2*d);
				surfaceArea = 0.125*d2*d2;
			}
		}
		
		// Now contact area is Vcell*S(d)/hperp
		surfaceArea *= Vcell/dist.x;
	}
	
	// times position if axisym
	if(fmobj->IsAxisymmetric()) surfaceArea *= ndptr->x;
	
	// output results
	*hperp = dist.x;
	return surfaceArea;
}

// retrieve -2*scale*(mass gradient) for material matfld
// and but set components zero on symmetry planes
// Do not call if volume gradient not used (or if volumeGradientIndex<0)
void CrackVelocityFieldMulti::GetVolumeGradient(int matfld,const NodalPoint *ndptr,Vector *grad,double scale) const
{
	CopyScaleVector(grad,&mvf[matfld]->contactInfo->terms[mpmgrid.volumeGradientIndex],scale);
	if(ndptr->fixedDirection&XSYMMETRYPLANE_DIRECTION) grad->x = 0.;
	if(ndptr->fixedDirection&YSYMMETRYPLANE_DIRECTION) grad->y = 0.;
	if(fmobj->IsAxisymmetric())
	{	// Need special case here to insure grad->z is zero (it is non-zero in volume gradient due to use of
		//   z component in extra axisymmetric shape function
		grad->z = 0.;
	}
	else if(ndptr->fixedDirection&ZSYMMETRYPLANE_DIRECTION)
    {	grad->z = 0.;
	}
}

#pragma mark VELOCITY METHODS

// Various calculations of grid values
// But do not called if it is a mirrored field (calculation in source will finish handle it)
void CrackVelocityFieldMulti::GridValueCalculation(int calcOption)
{	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidSourceField(mvf[i],fieldNum))
			mvf[i]->GridValueCalculation(calcOption);
	}
}

#pragma mark BOUNDARY CONDITIONS

// zero one component of moment and velocity on each material field
void CrackVelocityFieldMulti::ZeroVelocityBC(Vector *norm,int passType,double deltime,Vector *freaction)
{	for(int i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
        {	mvf[i]->ZeroVelocityBC(norm,passType,deltime,freaction);
        }
    }
}

// add one component momentum and velocity from BCs to each material field
void CrackVelocityFieldMulti::AddVelocityBC(Vector *norm,double vel,int passType,double deltime,Vector *freaction)
{	for(int i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
        {	mvf[i]->AddVelocityBC(norm,vel,passType,deltime,freaction);
        }
    }
}

// Reflect one component of velocity and momentum from a node
// This node is the BC node and rcvf is crack velocity field on the reflected node
// Could be issue is the two nodes have different arrangement of crack velocity fields
void CrackVelocityFieldMulti::ReflectVelocityBC(Vector *norm,CrackVelocityField *rcvf,double vel0,double reflectRatio,
												int passType,double deltime,Vector *freaction)
{	MatVelocityField **rmvf = rcvf->GetMaterialVelocityFields();
    for(int i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	if(MatVelocityField::ActiveNonrigidField(rmvf[i]))
			{	double rvel;
				if(passType==UPDATE_GRID_STRAINS_CALL)
					rvel = vel0 + reflectRatio*(vel0 - DotVectors(norm,rmvf[i]->vk));
				else
					rvel = vel0 + reflectRatio*(vel0 - DotVectors(norm,&rmvf[i]->pk)/rmvf[i]->mass);
				mvf[i]->AddVelocityBC(norm,rvel,passType,deltime,freaction);
			}
		}
	}
}

#pragma mark CONTACT ACCESSORS

// get total contact volume for all nonrigid materials
// if requireCracks is true than only count materials that see cracks, otherwise use all fields
double CrackVelocityFieldMulti::GetContactVolumeNonrigid(bool requireCracks) const
{	double volume = 0.;
	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidSeesCracksField(mvf[i],requireCracks))
		volume += mvf[i]->GetContactVolume();
	}
	return volume;
}

// get center of mass displacement (actually sum of displacement*mass so displacement is vector/total mass)
// Includes only non-rigid materials
// if requireCracks is true than only count materials that see cracks, otherwise use all fields
//     it is false in material contact, but true for crack contact
// If on symmetry plane, that component will be zeroed out
// Must have extrapolated term that is requested
Vector CrackVelocityFieldMulti::GetCMDisplacement(NodalPoint *np,bool requireCracks,bool useDisps) const
{	Vector dk;
	ZeroVector(&dk);
	if(useDisps)
	{	for(int i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveNonrigidSeesCracksField(mvf[i],requireCracks))
				AddVector(&dk,&mvf[i]->contactInfo->terms[mpmgrid.displacementIndex]);
		}
	}
	else
	{	for(int i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveNonrigidSeesCracksField(mvf[i],requireCracks))
				AddVector(&dk,&mvf[i]->contactInfo->terms[mpmgrid.positionIndex]);
		}
	}
	AdjustForSymmetry(np,&dk,false);
	return dk;
}

// get sum of forces for all material fields in this crack velocity field
Vector CrackVelocityFieldMulti::GetCMatFtot(void)
{	Vector fk;
	ZeroVector(&fk);
	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			AddVector(&fk,mvf[i]->GetFtotPtr());
	}
	return fk;
}

#pragma mark ACCESSORS

// Return true if this material field is physical field in memory
// Always true for NairnMPM
bool CrackVelocityFieldMulti::MVFInMemory(int mi) const
{	if(fieldNum==0 || !mvf[mi]->IgnoresCracks()) return true;
	// if here, the material velocity field is link to crack field [0] field
	return false;
}

// Look for presence on non rigid points
// This counts fields that do not see cracks, will need revise to look for nonrigid particles that see cracks
bool CrackVelocityFieldMulti::HasPointsNonrigid(void) const { return (numberPoints-numberRigidPoints) > 0; }

// Sum mass (only of non-rigid materials) and return the result
// If requireCracks is true then count only those that see cracks, otherwise
//		get entire mass
double CrackVelocityFieldMulti::GetTotalMass(bool requireCracks) const
{	double mass = 0;
	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidSeesCracksField(mvf[i],requireCracks))
			mass += mvf[i]->mass;
	}
	return mass;
}

// Count number of materials in each crack velocity field on this node
// Sum mass (only of non-rigid materials) and return the result
// Note: for mass, this does not count materials that ignore cracks unless field [0].
// Copy momenta to be restored after force extrapolation
// if has more than 1 material, set hasMaterialContact to true
// Uses: Only called once per times step in post extrapolation phase of mass an momentum task
double CrackVelocityFieldMulti::GetTotalMassAndCount(bool &hasMaterialContact)
{	
	double mass = 0.;
	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	numberMaterials++;
			if(!mvf[i]->IsRigidField())
			{   if(MVFInMemory(i))	mass += mvf[i]->mass;
				
				// copy the extrapolated momenta and contact XPIC if needed
				mvf[i]->GetTotalMassAndCount();
			}
		}
	}
	if(numberMaterials>1) hasMaterialContact = true;
	return mass;
}

// Return 1 or 0 if has any points that see cracks
int CrackVelocityFieldMulti::HasPointsThatSeeCracks(void)
{	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidSeesCracksField(mvf[i],true))
			return 1;
	}
	return 0;
}

// Return true or false if this crack velocity field has any field that ignores cracks
// This is only called for cvf field [0] on a node
bool CrackVelocityFieldMulti::HasFieldsThatIgnoreCracks(void) const
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	if(mvf[i]->IgnoresCracks())
				return true;
		}
	}
	return false;
}

// total mass and kinetic energy all velocity fields (rigid particles not counted)
// in g-mm^2/sec^s = nanoJ
void CrackVelocityFieldMulti::AddKineticEnergyAndMass(double &kineticEnergy,double &totalMass)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	if(mvf[i]->mass > 0.)
			{	totalMass += mvf[i]->mass;
				double magp = DotVectors(&mvf[i]->pk,&mvf[i]->pk);
				kineticEnergy += 0.5*magp/mvf[i]->mass;
			}
		}
	}
}

// total number of materials (rigid + nonrigid)
int CrackVelocityFieldMulti::GetNumberMaterials(void) { return numberMaterials; }

// Count the number of nonrigid materials
int CrackVelocityFieldMulti::GetNumberNonrigidMaterials(void)
{	int numMats = 0;
	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			numMats++;
	}
	return numMats;
}

// get center of mass momentum for all nonrigid material fields in this crack velocity field
//    or sum(pi) = sum(mi*vi*Sip) this crack velocity field
// Only counts materials that account for cracks
// If totalFtot (pre-zeroed) not NULL, get sum(fi) = sum(mi*ai*Sip) this crack velocity field
Vector CrackVelocityFieldMulti::GetCMatMomentum(bool &hasParticles,double *foundMass,Vector *totalFtot,bool useVelocity) const
{	Vector pk;
	ZeroVector(&pk);
	hasParticles = false;
	*foundMass = 0;
	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidSeesCracksField(mvf[i],true))
		{	if(useVelocity)
				AddScaledVector(&pk,&mvf[i]->vk[0],mvf[i]->mass);
			else
				AddVector(&pk,&mvf[i]->pk);
			*foundMass += mvf[i]->mass;
			hasParticles = true;
			if(totalFtot!=NULL) AddVector(totalFtot,mvf[i]->GetFtotPtr());
		}
	}
	return pk;
}

// add contact force on all rigid materials to the input vector
void CrackVelocityFieldMulti::SumAndClearRigidContactForces(Vector *fcontact,bool clearForces,double scale,Vector *ftotal)
{
	// if none, nothing to do
	if(numberRigidPoints==0) return;
	
	// check for rigid materials, but add only nonmirrored fields
	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveRigidField(mvf[i]))
		{   if(MVFInMemory(i))
			{	AddScaledVector(&fcontact[i],mvf[i]->GetFtotPtr(),scale);
				if(ftotal!=NULL) AddScaledVector(ftotal,mvf[i]->GetFtotPtr(),scale);
				if(clearForces) ZeroVector(mvf[i]->GetFtotPtr());
			}
		}
	}
}

/* in response to crack contact, change momentum by changing velocity of all 
	nonrigid materials the same amount (and only materials that account
	for cracks)
 
   Change velocity by dP/M, where M is total mass
   Material i velocity becomes vi = pi/mi + dP/M
   Material i momentum change is mi vi = pi + mi dP/M
*/
void CrackVelocityFieldMulti::ChangeCrackMomentum(Vector *delP,int callType,double deltime)
{
	int i;
	
#ifdef CHECK_NAN
	if(delP->x!=delP->x || delP->y!=delP->y || delP->z!=delP->z)
	{
#pragma omp critical (output)
		{	cout << "\n# CrackVelocityFieldMulti::ChangeCrackMomentum: bad delP";
			PrintVector(" = ",delP);
			cout << endl;
			Describe();
		}
	}
#endif
	
	// special case for only one material (and it must be nonrigid)
	if(numberMaterials==1)
	{	for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveField(mvf[i]))
			{	mvf[i]->ChangeMatMomentum(delP,callType,deltime);
				break;
			}
		}
	}
	
	// more than one material
	else
	{	Vector partialDelP;
		double totMass = GetTotalMass(true);
		for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveNonrigidSeesCracksField(mvf[i],true))
				mvf[i]->ChangeMatMomentum(CopyScaleVector(&partialDelP, delP, mvf[i]->mass/totMass),callType,deltime);
		}
	}
}

// copy all material velocity fields for boundary conditions methods, returning new offset into the save array
int CrackVelocityFieldMulti::CopyFieldMomenta(Vector *holdPk,int offset)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	holdPk[offset]=mvf[i]->pk;
			offset++;
		}
	}
	return offset;
}

#if ADJUST_COPIED_PK == 1
// set symmetry plane momenta to zero
void CrackVelocityFieldMulti::AdjustForSymmetryBC(NodalPoint *ndptr)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			mvf[i]->AdjustForSymmetryBC(ndptr->fixedDirection);
	}
}
#endif

// paste all material velocity fields back for boundary conditions methods, returning new offset into the saved array
int CrackVelocityFieldMulti::PasteFieldMomenta(Vector *holdPk,int offset)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	mvf[i]->pk=holdPk[offset];
			offset++;
		}
	}
	return offset;
}

// for debugging
void CrackVelocityFieldMulti::Describe(void) const
{	CrackVelocityField::Describe();
	cout << "#     multimaterial: nmat= " << numberMaterials << " nrigidpts= " << numberRigidPoints << endl;
	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			mvf[i]->Describe(i);
	}
}
