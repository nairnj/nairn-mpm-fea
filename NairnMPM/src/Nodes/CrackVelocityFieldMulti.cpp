/********************************************************************************
	CrackVelocityFieldMulti.cpp
	nairn-mpm-fea

	Created by John Nairn on 21 August 2009.
	Copyright (c) 2009 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Nodes/CrackVelocityFieldMulti.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Materials/ContactLaw.hpp"
#include "Nodes/MaterialContactNode.hpp"
#include "Exceptions/MPMWarnings.hpp"

#define MASS_MIN 1.e-5
#define MASS_MAX 0.99999

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
		{   // Always delete
			delete mvf[i];
		}
	}
}

// zero all active material velocity fields
// but mirrored fields just set to NULL
void CrackVelocityFieldMulti::ZeroMatFields(void)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(mvf[i]!=NULL)
		{	if(mvf[i]->numberPoints>0)
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

// Add to mass gradient
// This gradient is only used in multimaterial contact calculations
// Corrections for symmetry planes and axisymmetry are done in GetVolumeGradient(). It is
//		faster there because that is called less while this is called for every node-particle pair
void CrackVelocityFieldMulti::AddVolumeGradient(int matfld,MPMBase *mptr,double dNdx,double dNdy,double dNdz)
{	double Vp = mptr->GetVolume(DEFORMED_AREA_FOR_GRADIENT);
	mvf[matfld]->volumeGrad->x += Vp*dNdx;
	mvf[matfld]->volumeGrad->y += Vp*dNdy;
	mvf[matfld]->volumeGrad->z += Vp*dNdz;
}

// call this when copy from ghost to real node to sume gradients
void CrackVelocityFieldMulti::CopyVolumeGradient(int matfld,Vector *grad)
{	mvf[matfld]->volumeGrad->x += grad->x;
	mvf[matfld]->volumeGrad->y += grad->y;
	mvf[matfld]->volumeGrad->z += grad->z;
}

// This is called for crack field [i]. Each material field in crack source [0] that ignores cracks should be copied
//    to the corresponding field in this crack velocity field
// For NairnMPM, a new field is created if needed (instead of just mirroring it) and then filled with copy
void CrackVelocityFieldMulti::MirrorFieldsThatIgnoreCracks(MatVelocityField *rmvf,int rigidFieldNum)
{	
	// create material field in this crack velocity field if needed, otherwise, just be sure it is zeroed.
	if(mvf[rigidFieldNum]==NULL)
	{	mvf[rigidFieldNum] = new MatVelocityField(RIGID_FIELD_BIT);
		if(mvf[rigidFieldNum]==NULL) throw CommonException("Memory error allocating material velocity field.",
													"CrackVelocityFieldMulti::MirrorFieldsThatIgnoreCracks");
	}
	else
		mvf[rigidFieldNum]->Zero();
 	
	// add number rigid points this crack velocity field
	numberRigidPoints += rmvf->numberPoints;
	numberPoints += rmvf->numberPoints;
	
	// copy all extrapolated items
    
    // momentum, number of points, and velocity
	Vector rvel = rmvf->GetVelocity();
    mvf[rigidFieldNum]->AddMomentumTask1(&rmvf->pk,&rvel,rmvf->numberPoints);
    
    // mass and volume
	mvf[rigidFieldNum]->mass = rmvf->mass;
	mvf[rigidFieldNum]->AddContactVolume(rmvf->GetContactVolume());
    
    // displacement and volume gradient
	CopyVector(&mvf[rigidFieldNum]->disp,&rmvf->disp);
	CopyVector(mvf[rigidFieldNum]->volumeGrad,rmvf->volumeGrad);
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
            {	ZeroVector(&mvf[i]->pk);
                ZeroVector(&mvf[i]->disp);
                if(mvf[i]->volumeGrad!=NULL) ZeroVector(mvf[i]->volumeGrad);
				mvf[i]->SetContactVolume(0.);
            }
            else
            {   // for rigid particles, keep initial pk
                // can project displacement using current velocity because
                // particle mass is its volume
                // dnew = Sum (Vp*fpi*(d + v dt)) = dold + Sum (Vp*fpi*v*dt) = dold + pk*dt
                AddScaledVector(&mvf[i]->disp,&mvf[i]->pk,deltaTime);
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
void CrackVelocityFieldMulti::UpdateMomentaOnField(double timestep)
{	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidSourceField(mvf[i],fieldNum))
			mvf[i]->UpdateMomentum(timestep);
    }
}

#pragma mark MATERIAL CONTACT

/* Called in multimaterial mode to check contact at nodes with multiple materials

	Input parameters:
		mvf[]->mass,pk,volumeGrad,disp (if contact by displacements)

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
		RigidMaterialContactOnCVF(rigidMat,multiRigid,mcn,deltime,callType);
	else
	{
		MaterialContactOnCVFLumped(mcn,deltime,callType,activeMat,numMats,Pc,Mc);
	}
	
	delete [] activeMat;
	return;
}

// Contact for each material lumped with the other materials
void CrackVelocityFieldMulti::MaterialContactOnCVFLumped(MaterialContactNode *mcn,double deltime,int callType,int *activeMat,int numMats,
														 Vector Pc,double Mc)
{
	// get the node
	NodalPoint *ndptr = mcn->GetTheNode();
	
	double qrate = 0.;
	int i;
	
	// loop materials or pairs of materials
	bool doingPairs = (numMats==2 && contact.materialNormalMethod!=EACH_MATERIALS_MASS_GRADIENT);	// has only one pair;
	int miMax = doingPairs ? numMats-1 : numMats ;
	
	// Get center of mass displacememnt used in this mode
	Vector dispc = GetCMDisplacement(ndptr,false);
	ScaleVector(&dispc,1./Mc);
	AdjustForSymmetry(ndptr,&dispc,false);
	
	MatVelocityField *mvfi,*mvfj;
	for(int mi=0;mi<miMax;mi++)
	{
		// known active field
		i = activeMat[mi];
		mvfi = mvf[i];
		double massi = mvfi->mass;
		double voli = mvfi->GetContactVolume();
		double volj = GetVolumeNonrigid(false) - voli;
		
		// Contact law from other material with most volume
		// Also find volume gradient of all other non-rigid materials
		int j = -1;
		Vector gradj;
		ZeroVector(&gradj);
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
			Vector normj;
			GetVolumeGradient(jj,ndptr,&normj,-1.);
			AddVector(&gradj,&normj);
		}
		
		// problem if paired material not found, but it will always be found
		if(j<0) continue;
		mvfj = mvf[j];
		ContactLaw *theContactLaw = contact.GetMaterialContactLaw(i,j);
		
		// compared to lumped mass (or other material)
		double mred = (Mc-massi)/Mc;				// final mred = massi*(Mc-massi)/Mc, but found later when needed
		double massRatio = massi/Mc;
		
		// some variables
		Vector norm,delPi,delta,tangDel;
		double dotn,deln=0.,delt=0.;				// normal force and delta components
		double contactArea = -1.,mredDelWf = -1.,contactGridN = 1.;
		bool hasDeln;
			
		// Determine if in contact
		bool comContact = false;
		
		// find mi(vc-vi) = (mi/mc)pc-pi = (mi pc - mc pi)/mc or momentum change to match ctr of mass momentum
		CopyScaleVector(&delPi,&mvfi->pk,-1.);
		AddScaledVector(&delPi,&Pc,massRatio);
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

			// Get normal vector, dotn, and deln
			// Also find displacment vector
			
			// Get normal vector by various options
			hasDeln = false;
			if(!NonRigidCustomNormal(ndptr,i,j,norm))
			{	norm = GetNormalVector(mcn,i,-1,voli,&gradj,volj,1.,hasDeln,deln);
				
				// nan here means volume gradients zero as gradient normal to a symmetry plane
				if(norm.x!=norm.x || norm.y!=norm.y || norm.z!=norm.z) continue;
			}
			
			// Get sticking force (actual (vc-vi) = delPi/mi and (vb-v1) = delPi/mred (when mred is done)
			dotn=DotVectors(&delPi,&norm);
			
			// Get COD and displacement of each material and final reduced mass
			// Let dispc = xiCOM, dispa = xia = mvfi->disp/massi, and mred = (Mc-Mi)/Mc
			// From my 2013 paper, COD = xib - xia = (Mc/(Mc-mi))(xiCOM - xia) = (dispc-xia)/mred
			Vector dispa;
			CopyScaleVector(&dispa,&mvfi->disp,1./massi);
			AdjustForSymmetry(ndptr,&dispa,false);
			CopyScaleVector(&delta,&dispc,1./mred);
			AddScaledVector(&delta,&dispa,-1./mred);
			
			// Final reduced mass = (mi*(Mc-mi))/Mc while until here was (Mc-mi)/Mc
			mred *= massi;
			
			// get deln = delta.n corrected if simulation is extrapolating positions or displacements
			if(hasDeln)
			{	// accept tangential component of delta, but change normal to deln
				double extrapDeln = DotVectors(&delta,&norm);
				AddScaledVector(&delta,&norm,deln-extrapDeln);
			}
			else
			{	// get dispb = delta + xia = xib-xia+xia = xib (already adjusted for symmetry)
				Vector dispb = delta;
				AddVector(&dispb,&dispa);
				deln = contact.MaterialSeparation(&delta,&dispa,&dispb,&norm,ndptr);
			}
			
		}
			
		// Here inContact is usually true, but may be false for imperfet interface or for
		// any other contact law that want to continue (e.g., friction with adhesion)
		if(comContact)
		{	// Use delPi found above and no further change needed, no other data needed either
		}
		else if(theContactLaw->IsFrictionalContact())
		{	// Handle frictional law
			
			// get area if needed this friction law
			if(contactArea<0. && theContactLaw->ContactLawNeedsContactArea())
			{	contactArea = GetContactArea(ndptr,voli,volj,
											 &delta,&norm,&tangDel,&delt,&contactGridN);
			}
			
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
														  contactArea,deltime,delFiPtr))

			{	continue;
			}
		}
		else
		{	// Handle imperfect interface
			
			// Interfaces always need contact area and tangDel and delt (find if not done already)
			if(contactArea<0.)
			{	contactArea = GetContactArea(ndptr,voli,volj,
											 &delta,&norm,&tangDel,&delt,&contactGridN);
			}
			
			// Find delFi = (ma Fb - mb Fa)/Mc = (ma Fc - Mc Fa)/Mc = (ma Fc/Mc) - Fa (when needed)
			Vector fImp;
			if(callType==UPDATE_MOMENTUM_CALL)
			{	fImp = GetCMatFtot();
				ScaleVector(&fImp, massi/Mc);
				AddScaledVector(&fImp, mvfi->GetFtotPtr(), -1.);
			}
			else
				ZeroVector(&fImp);
			
			// Get interface force and energy
			double rawEnergy;
			theContactLaw->GetInterfaceForces(&norm,&fImp,&rawEnergy,contactArea,&delPi,dotn,mred,
											  &tangDel,deln,delt,contactGridN);
			if(callType==UPDATE_MOMENTUM_CALL)
			{	// add force (if any) to momentum change
				AddScaledVector(&delPi, &fImp, timestep);
				
				// Add interface energy. (Legacy units g-mm^2/sec^2 or multiply by 1e-9 to get J - kg-m^2/sec^2)
				NodalPoint::interfaceEnergy += rawEnergy;
			}
		}
		
		// ----------------------------------------------------------------------------------
		// 7. Apply momentum change to the velocity fields
		
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
void CrackVelocityFieldMulti::RigidMaterialContactOnCVF(int rigidFld,bool multiRigid,MaterialContactNode *mcn,double deltime,int callType)
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
	{	if(!MatVelocityField::ActiveField(mvf[i]) || mvf[i]->IsRigidField()) continue;
		
		// First determine contact law with rigid material
		ContactLaw *theContactLaw = contact.GetMaterialContactLaw(i,rigidFld);
		
		// get contact parameters
		Vector delta,tangDel;
		Vector dispRigid,dispi;				// displacement of rigid material and this material
		double deln=0.,delt=0.;				// delta components
		bool hasDeln;
 	
		// some variables
		Vector norm,delPi;
		double dotn,massi=mvf[i]->mass,voli=mvf[i]->GetContactVolume();
		double contactArea = -1.,mredDelWf = -1.,contactGridN = 1.;
		
		// find mi(vr-vi) = - pi + mi*vr or momentum change to match rigid particle velocity
		CopyScaleVector(&delPi,&mvf[i]->pk,-1.);
		AddScaledVector(&delPi,&rvel,massi);
		AdjustForSymmetry(ndptr,&delPi,false);
		
		// if needed, determine if the surfaces are in contact
		bool comContact = false;
		if(theContactLaw->IgnoreContact() || (fmobj->dflag[6]&4))
		{   // will use delPi to match rigid particle velocity
			comContact = true;
		}
		
		else
		{   // Ignore very small interactions
			double volRatio=voli/(rigidVolume+GetVolumeNonrigid(false));
			if(volRatio<MASS_MIN || volRatio>MASS_MAX)
			{	continue;
			}
			
			// Get normal vector, dotn, and deln
			// Also find displacment vector
			
			// Get normal vector by various options
			hasDeln = false;
			if(!RigidCustomNormal(ndptr,i,rigidFld,norm))
			{	Vector gradj;
				GetVolumeGradient(rigidFld,ndptr,&gradj,-1.);
				norm = GetNormalVector(mcn,i,rigidFld,voli,&gradj,rigidVolume,contact.rigidGradientBias,hasDeln,deln);
				
				// nan here means volume gradients zero as gradient normal to a symmetry plane
				if(norm.x!=norm.x || norm.y!=norm.y || norm.z!=norm.z) continue;
			}
			
			// get approach direction momentum from delPi.n (actual (vc-vi) = delPi/mi)
			dotn=DotVectors(&delPi,&norm);
			
			// Displacement calculations: rigid material displacement was scaled by volume, while non-rigid was weighted by mass
			CopyScaleVector(&dispRigid,&mvf[rigidFld]->disp,1./actualRigidVolume);
			AdjustForSymmetry(ndptr,&dispRigid,false);
			CopyScaleVector(&dispi,&mvf[i]->disp,1./massi);
			AdjustForSymmetry(ndptr,&dispi,false);
			delta = dispRigid;
			SubVector(&delta,&dispi);
			
			// get deln = delta.n corrected if simulation is extrapolating positions or displacements
			if(hasDeln)
			{	// accept tangential component of delta, but change normal to deln
				double extrapDeln = DotVectors(&delta,&norm);
				AddScaledVector(&delta,&norm,deln-extrapDeln);
			}
			else
			{	deln = contact.MaterialSeparation(&delta,&dispi,&dispRigid,&norm,ndptr);
			}
		}
		
		// Here inContact is usually true, but may be false for imperfect interface or for
		// any other contact law that wants to continue (e.g., friction with adhesion)
		if(comContact)
		{	// Use delPi found above and no further change needed, no other data needed either
		}
		else if(theContactLaw->IsFrictionalContact())
		{	// Handle frictional law
			
			// may need area and contactGridN, but get contactGridN even if do not need area
			if(contactArea<0. && theContactLaw->ContactLawNeedsContactArea())
			{	contactArea = GetContactArea(ndptr,voli,rigidVolume,&delta,&norm,&tangDel,&delt,&contactGridN);
			}
			else
			{	Vector gridDist = mpmgrid.GetPerpendicularDistance(&norm,ndptr);
				contactGridN = gridDist.x;
			}
			
			// Find dF = (ma Fb - mb Fa)/Mc = -Fa for rigid material (infinite mb)
			// note that dacc = dF/mred = -Fa/ma
			bool getHeating = false;
			Vector delFi;
			Vector *delFiPtr = NULL;
			if(callType==UPDATE_MOMENTUM_CALL)
			{	delFi = mvf[i]->GetFtot();
				ScaleVector(&delFi, -1.);
				delFiPtr = &delFi;
				// second order heating needs acceleration too
				if(ConductionTask::matContactHeating) getHeating = true;
			}
			
			// get change in momentum, but false return means friction law wants no change (i.e, continue as if not in contact)
			if(!theContactLaw->GetFrictionalDeltaMomentum(&delPi,&norm,dotn,deln,&mredDelWf,massi,getHeating,
														  contactArea,deltime,delFiPtr))
			{	continue;
			}
		}
		else
		{	// Handle imperfect interface
			
			// Interfaces always need contact area (find if not done already)
			if(contactArea<0.)
			{	contactArea = GetContactArea(ndptr,voli,rigidVolume,&delta,&norm,&tangDel,&delt,&contactGridN);
			}
			
			// Find delFi = (ma Fb - mb Fa)/Mc = (ma Fc - Mc Fa)/Mc = (ma Fc/Mc) - Fa = -Fa (when needed)
			Vector fImp;
			ZeroVector(&fImp);
			if(callType==UPDATE_MOMENTUM_CALL)
				AddScaledVector(&fImp, mvf[i]->GetFtotPtr(), -1.);
			
			// Get interface force and energy
			double rawEnergy;
			theContactLaw->GetInterfaceForces(&norm,&fImp,&rawEnergy,
													contactArea,&delPi,dotn,massi,&tangDel,deln,delt,contactGridN);
			if(callType==UPDATE_MOMENTUM_CALL)
			{	// add force (if any) to momentum change
				AddScaledVector(&delPi, &fImp, timestep);
				
				// Add interface energy. (Legacy units g-mm^2/sec^2 or multiply by 1e-9 to get J - kg-m^2/sec^2)
				NodalPoint::interfaceEnergy += rawEnergy;
			}
		}
		
		// ----------------------------------------------------------------------------------
		// 6. Apply momentum change to the velocity fields
		
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
				double qrate = mredDelWf/massi;
				conduction->AddFluxCondition(ndptr,fabs(qrate/deltime),true);
				
#pragma omp atomic
				NodalPoint::frictionWork += qrate;
			}
		}
	}
}

// Get normal vector for material i to j by various options
// Note 1: j<0 means others lumped (info passed to LR methods, other methods don't used j)
// Note 2: gradJ and otherVolume may be undefined if using LR method (because not used)
// Note 3: If LR method, hasDeln with change to true and deln will be set to normal separation
Vector CrackVelocityFieldMulti::GetNormalVector(MaterialContactNode *mcn,int i,int j,double voli,Vector *gradj,double otherVolume,
												double jBias,bool &hasDeln,double &deln)
{
	NodalPoint *ndptr = mcn->GetTheNode();
	Vector norm;
	
	// Pick and option
	switch(contact.materialNormalMethod)
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
			norm = contact.contactNormal;
			AdjustForSymmetry(ndptr,&norm,true);
			break;
			
		default:
			break;
	}
	
	// final result
	return norm;
}

// For contact between nonrigid material i and j, see if any dflag[] normals are
// current in use
bool CrackVelocityFieldMulti::NonRigidCustomNormal(NodalPoint *ndptr,int i,int j,Vector &norm)
{
	//----------------------------------------------------------------------------------
	// Various developer flag hacks to change the normal vector in specific problems
	// current options: 5 (spherical particle, only here in nonrigid)
	if(fmobj->dflag[0]==5)
	{	// Radial normal for spherical particle at the origin
		norm.x = ndptr->x;
		norm.y = ndptr->y;
		norm.z = ndptr->z;
		double normmag=sqrt(DotVectors(&norm,&norm));
		norm.x /= normmag;
		norm.y /= normmag;
		norm.z /= normmag;
		AdjustForSymmetry(ndptr,&norm,true);
		return true;
	}
	else if(fmobj->dflag[0]==4 && (MaterialBase::GetFieldMatID(i)==0 || MaterialBase::GetFieldMatID(j)==0))
	{	// NOTE: if called for rigid contact, that code will replace these calculations with new normal
		// use special normals for cutting simulation with rake angle in dflag[1]
		// and the material below the crack as the first defined material
		// Assumes cut material is 1 (ID=0) and tool material is 2 (ID 1)
		if(MaterialBase::GetFieldMatID(i)==1 || MaterialBase::GetFieldMatID(j)==1)
		{	Vector nrpos;
			// position of material being cut
			if(MaterialBase::GetFieldMatID(i)==0)
				CopyScaleVector(&nrpos,&mvf[i]->disp,1./mvf[i]->mass);
			else
				CopyScaleVector(&nrpos,&mvf[j]->disp,1./mvf[j]->mass);
			if(fmobj->dflag[1]>-90.)
			{   // use rake angle as specified in dflag[1]
				if(nrpos.y<0.)
				{	// from below the cut line up to tool
					norm.x=0.;
					norm.y=1.;
				}
				else
				{	double radAngle=(double)fmobj->dflag[1]*PI_CONSTANT/180.;
					norm.x=cos(radAngle);
					norm.y=-sin(radAngle);
				}
				// Flip if needed
				if(MaterialBase::GetFieldMatID(i)==1)
				{	norm.x = -norm.x;
					norm.y = -norm.y;
				}
			}
			else
			{	// When rake angle < -90, use the normal calculated by standard methods
				// An option in cutting is to apply normal on rubbing surface
				//		as vertical (with following code) or to continue
				//		with the caculated one (if commented out)
				/*
				 if(nrpos.y<0.)
				 {	// from below the cut line up to tool
				 norm.x=0.;
				 norm.y = MaterialBase::GetFieldMatID(i)==1 ? -1. : 0.;
				 }
				 */
			}
		}
		
		norm.z=0.;
		AdjustForSymmetry(ndptr,&norm,true);
		return true;
	}
	
	// not custom normal available
	return false;
}

// For contact between nonrigid material i and rigidFld, see if any dflag[] normals are
// currently in use
bool CrackVelocityFieldMulti::RigidCustomNormal(NodalPoint *ndptr,int i,int rigidFld,Vector &norm)
{
	// Various developer flag hacks to change the normal vector in specific problems
	// current options: 4 (cutting, but only here in rigid contact)
	if(fmobj->dflag[0]==4 && MaterialBase::GetFieldMatID(i)==0)
	{	if(MaterialBase::GetFieldMatID(rigidFld)==1)
		{	// use special normals for cutting simulation with rake angle in dflag[1]
			// and the material below the crack as the first defined material
			// Here material 1 (ID=0) = material to cut and  2 (ID=1) is tool
			Vector nrpos;
			CopyScaleVector(&nrpos,&mvf[i]->disp,1./mvf[i]->mass);
			if(fmobj->dflag[1]>-90.)
			{   // use rake angle as specified in dflag[1]
				if(nrpos.y<0.)
				{	// from below the cut line up to tool
					norm.x=0.;
					norm.y=1.;
				}
				else
				{	double radAngle=(double)fmobj->dflag[1]*PI_CONSTANT/180.;
					norm.x=cos(radAngle);
					norm.y=-sin(radAngle);
				}
				norm.z=0.;
				AdjustForSymmetry(ndptr,&norm,true);
			}
			else
			{	// When rake angle < -90, use the normal calculated by standard methods
				// An option in cutting is to apply normal on rubbing surface
				//		as vertical (with following code) or to continue
				//		with the caculated one (if commented out)
				if(nrpos.y<0.)
				{	// from below the cut line up to tool
					norm.x=0.;
					norm.y=1.;
				}
				norm.z=0.;
				AdjustForSymmetry(ndptr,&norm,true);
			}
			
		}
		else if(MaterialBase::GetFieldMatID(rigidFld)==2 && fmobj->dflag[2]>0)
		{	// Here material 1 (ID=0) = material to cut and material 3 (ID=2) is plane base
			// must specify thickness in dflag[2] in microns
			Vector nrpos;
			CopyScaleVector(&nrpos,&mvf[i]->disp,1./mvf[i]->mass);
			norm.x = 0.;
			int num = ndptr->num;
			double celly = nd[num+mpmgrid.yplane]->y-ndptr->y;
			norm.y = nrpos.y < ((double)fmobj->dflag[2]/1000.)+1.5*celly ? 1. : -1. ;			// from below or above the plane
			norm.z=0.;
			AdjustForSymmetry(ndptr,&norm,true);
		}
		
		return true;
	}
	
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

// Get interfacial contact area and also find tangDel and delt (only used in interface laws)
// Input is ndoptr, voli, volb, delta, norm
// Ouput is tangDel, delt, hperp
double CrackVelocityFieldMulti::GetContactArea(NodalPoint *ndptr,double voli,double volb,Vector *delta,
											   Vector *norm,Vector *tangDel,double *deltout,double *hperp) const
{
    // tangential vector in tangDel (might be zero vector if delta || to n)
	double deln = DotVectors(delta,norm);
    CopyVector(tangDel,delta);
    AddScaledVector(tangDel,norm,-deln);				// delta - deln (n) = delt (t)
	
	// by normalizing to positive delt, hat t always points in positive direction
    double delt=sqrt(DotVectors(tangDel,tangDel));
    if(!DbleEqual(delt,0.)) ScaleVector(tangDel,1./delt);
	
	// perpendicular distance to correct contact area and contact by positions
	Vector dist = mpmgrid.GetPerpendicularDistance(norm,ndptr);
	
	// Get raw surface area, it is divided by hperp to get contact area
	// Scale voltot=voli+volb to voltot*sqrt(2*vmin/voltot) = sqrt(2*vmin*vtot)
	// dist weightings to allow for Tartan grid
	// multiple by position for axisymmetric
	double surfaceArea = sqrt(2.*fmin(voli*dist.y,volb*dist.z)*(voli+volb))/dist.x;
	if(fmobj->IsAxisymmetric()) surfaceArea *= ndptr->x;			// times position if axisym
	
	// output components
	*deltout = delt;
	*hperp = dist.x;
	return surfaceArea;
}

// retrieve volume gradient for matnum (1 based) in crack field only (or zero if
// not there or not tracked
bool CrackVelocityFieldMulti::HasVolumeGradient(int matfld) const
{	if(mvf[matfld]==NULL) return false;
	if(mvf[matfld]->volumeGrad==NULL) return false;
	return true;
}

// retrieve -2*scale*(mass gradient) for material matfld
// and but set components zero on symmetry planes
void CrackVelocityFieldMulti::GetVolumeGradient(int matfld,const NodalPoint *ndptr,Vector *grad,double scale) const
{
	if(fmobj->IsAxisymmetric())
	{	// Need special case here to insure grad->z is zero (it is non-zero in volumeGrad due to use of
		//   z component in extra axisymmetric shape function
		if(ndptr->fixedDirection&XSYMMETRYPLANE_DIRECTION)
		{	grad->x=0.;
			grad->y = scale*mvf[matfld]->volumeGrad->y ;
		}
		else if(ndptr->fixedDirection&YSYMMETRYPLANE_DIRECTION)
		{	grad->y = 0.;
			grad->x = scale*mvf[matfld]->volumeGrad->x ;
		}
		else
		{	grad->x = scale*mvf[matfld]->volumeGrad->x;
			grad->y = scale*mvf[matfld]->volumeGrad->y;
		}
		grad->z=0.;
	}
	else
    {   CopyScaleVector(grad,mvf[matfld]->volumeGrad,scale);
		if(ndptr->fixedDirection&XSYMMETRYPLANE_DIRECTION) grad->x = 0.;
		if(ndptr->fixedDirection&YSYMMETRYPLANE_DIRECTION) grad->y = 0.;
		if(ndptr->fixedDirection&ZSYMMETRYPLANE_DIRECTION) grad->z = 0.;
	}
}

#pragma mark VELOCITY METHODS

// Calculate velocity at a node from current momentum and mass matrix in all velocity fields
void CrackVelocityFieldMulti::CalcVelocityForStrainUpdate(void)
{	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidSourceField(mvf[i],fieldNum))
			mvf[i]->CalcVelocityForStrainUpdate();
	}
}

#pragma mark BOUNDARY CONDITIONS

// zero one component of moment and velocity on each material field
void CrackVelocityFieldMulti::SetMomVel(Vector *norm,int passType)
{	for(int i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
        {	mvf[i]->SetMomentVelocityDirection(norm,passType);
        }
    }
}

// add one component momentum and velocity from BCs to each material field
void CrackVelocityFieldMulti::AddMomVel(Vector *norm,double vel,int passType)
{	for(int i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
        {	mvf[i]->AddMomentVelocityDirection(norm,vel,passType);
        }
    }
}

// Reflect one component of velocity and momentum from a node
// This node is the BC node and rcvf is crack velocity field on the reflected node
// Could be issue is the two nodes have different arrangement of crack velocity fields
void CrackVelocityFieldMulti::ReflectMomVel(Vector *norm,CrackVelocityField *rcvf,double vel0,double reflectRatio,int passType)
{	MatVelocityField **rmvf = rcvf->GetMaterialVelocityFields();
    for(int i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	if(MatVelocityField::ActiveNonrigidField(rmvf[i]))
			{	double rvel = vel0 + reflectRatio*(vel0 - DotVectors(norm,&rmvf[i]->pk)/rmvf[i]->mass);
				mvf[i]->AddMomentVelocityDirection(norm,rvel,passType);
			}
		}
	}
}

// set force in direction norm to -p(interpolated)/time such that updated momentum
//    of pk.i + deltime*ftot.i will be zero along norm
void CrackVelocityFieldMulti::SetFtotDirection(Vector *norm,double deltime,Vector *freaction)
{	for(int i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
        {	mvf[i]->SetFtotDirection(norm,deltime,freaction);
        }
    }
}

// add one component of force such that updated momentum will be mass*velocity
void CrackVelocityFieldMulti::AddFtotDirection(Vector *norm,double deltime,double vel,Vector *freaction)
{	for(int i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
        {	mvf[i]->AddFtotDirection(norm,deltime,vel,freaction);
        }
    }
}

// add one component of force such that updated momentum will be mass*velocity
void CrackVelocityFieldMulti::ReflectFtotDirection(Vector *norm,double deltime,CrackVelocityField *rcvf,
												   double vel0,double reflectRatio,Vector *freaction)
{	MatVelocityField **rmvf = rcvf->GetMaterialVelocityFields();
    for(int i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	if(MatVelocityField::ActiveNonrigidField(rmvf[i]))
			{	double rvel = vel0 + reflectRatio*(vel0 - DotVectors(norm,&rmvf[i]->pk)/rmvf[i]->mass);
				mvf[i]->AddFtotDirection(norm,deltime,rvel,freaction);
			}
		}
	}
}

#pragma mark ACCESSORS

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
// Uses: Only called once per times step in post extrapolation phase of mass an momentum task
// if has more than 1 material, set hasMaterialContact to true
double CrackVelocityFieldMulti::GetTotalMassAndCount(bool &hasMaterialContact)
{	
	double mass = 0.;
	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	numberMaterials++;
			if(!mvf[i]->IsRigidField())
			{	mass += mvf[i]->mass;
				
				// copy the extrapolated momenta
				mvf[i]->xpic[MatVelocityField::pkCopy] = mvf[i]->pk;
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

// get total contact volume for all nonrigid materials
// if requireCracks is true than only count materials that see cracks, otherwise use all fields
double CrackVelocityFieldMulti::GetVolumeNonrigid(bool requireCracks) const
{	double volume = 0.;
	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidSeesCracksField(mvf[i],requireCracks))
			volume += mvf[i]->GetContactVolume();
	}
	return volume;
}

// get total contact volume for all materials (including those that ignore cracks)
//     contact volume is actual volume except it is area in axisymmetry
// WARNING: this doubles the volume for symmetry planes (e.g. r=0 in axisymmetry
//   to enable volume screening to work
double CrackVelocityFieldMulti::GetVolumeTotal(NodalPoint *ndptr) const
{
	double volume = 0.;
	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			volume += mvf[i]->GetContactVolume();
	}
	
	// correct for on symmetry plane(s)
	if(ndptr!=NULL)
	{	if(ndptr->fixedDirection&XSYMMETRYPLANE_DIRECTION) volume *= 2.;
		if(ndptr->fixedDirection&YSYMMETRYPLANE_DIRECTION) volume *= 2.;
		if(ndptr->fixedDirection&ZSYMMETRYPLANE_DIRECTION) volume *= 2.;
	}
	return volume;
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
Vector CrackVelocityFieldMulti::GetCMatMomentum(bool &hasParticles,double *foundMass,Vector *totalFtot) const
{	Vector pk;
	ZeroVector(&pk);
	hasParticles = false;
	*foundMass = 0;
	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidSeesCracksField(mvf[i],true))
		{	AddVector(&pk,&mvf[i]->pk);
			*foundMass += mvf[i]->mass;
			hasParticles = true;
			if(totalFtot!=NULL) AddVector(totalFtot,mvf[i]->GetFtotPtr());
		}
	}
	return pk;
}

// get center of mass displacement (actually sum of displacement*mass so displacement is vector/total mass)
// Includes only non-rigid materials
// if requireCracks is true than only count materials that see cracks, otherwise use all fields
// If on symmetry plane, that component will be zeroed out
Vector CrackVelocityFieldMulti::GetCMDisplacement(NodalPoint *np,bool requireCracks) const
{	Vector dk;
	ZeroVector(&dk);
	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidSeesCracksField(mvf[i],requireCracks))
			AddVector(&dk,&mvf[i]->disp);
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

// add contact force on all rigid materials to the input vector
void CrackVelocityFieldMulti::SumAndClearRigidContactForces(Vector *fcontact,bool clearForces,double scale,Vector *ftotal)
{
	// if none, nothing to do
	if(numberRigidPoints==0) return;
	
	// check for rigid materials, but add only nonmirrored fields
	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveRigidField(mvf[i]))
		{	if(fieldNum==0)
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

#ifdef ADJUST_EXTRAPOLATED_PK_FOR_SYMMETRY
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

// get first active rigid field or return NULL. Also return number in rigidFieldNum
// Onluy called when copying rigid filed to field [0]
MatVelocityField *CrackVelocityFieldMulti::GetRigidMaterialField(int *rigidFieldNum)
{
	// if none return NULL
	if(numberRigidPoints==0) return NULL;
	
	// find the rigid field
	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveRigidField(mvf[i]))
		{	*rigidFieldNum=i;
			return mvf[i];
		}
	}
	
	return NULL;
}




