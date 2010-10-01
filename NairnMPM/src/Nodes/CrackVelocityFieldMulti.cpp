/********************************************************************************
	CrackVelocityFieldMulti.cpp
	NairnMPM

	Created by John Nairn on 21 August 2009.
	Copyright (c) 2009 John A. Nairn, All rights reserved.
********************************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Nodes/CrackVelocityFieldMulti.hpp"
#include "Exceptions/MPMTermination.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Materials/MaterialBase.hpp"
//#include "MPM_Classes/MPMBase.hpp"

#pragma mark INITIALIZATION

// constructor
CrackVelocityFieldMulti::CrackVelocityFieldMulti(short theLoc,int cnum) : CrackVelocityField(theLoc,cnum)
{	numberMaterials=0;
	numberRigidPoints=0;
}

// zero all active material velocity fiels
void CrackVelocityFieldMulti::ZeroMatFields(void)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			mvf[i]->Zero();
	}
	numberMaterials=0;
	numberRigidPoints=0;
}

#pragma mark TASK 1 METHODS

// add "mass" for rigid particle (task 1) - it counts particles, but mass will be zero
void CrackVelocityFieldMulti::AddMassTask1(int matfld) { numberRigidPoints++; }

// Add to mass gradient
void CrackVelocityFieldMulti::AddMassGradient(int matfld,double mp,double dNdx,double dNdy,double dNdz)
{	mvf[matfld]->massGrad->x+=mp*dNdx;
	mvf[matfld]->massGrad->y+=mp*dNdy;
	mvf[matfld]->massGrad->z+=mp*dNdz;
}

// Delete empty velocity fields, count number of materials, and return total mass
double CrackVelocityFieldMulti::GetTotalMassAndCount(void)
{	int i;
	double mass=0.;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	numberMaterials++;
			mass+=mvf[i]->mass;
		}
	}
	return mass;
}

// copy rigid material from another velocity field (cvfm) and add to mvf[rigidFieldNum] in this cvf
// This is only called if COMBINE_RIGID_MATERIALS is defined
void CrackVelocityFieldMulti::CombineRigidFrom(CrackVelocityFieldMulti *cvfm,int rigidFieldNum)
{
	// get other field, exit if none, or error if different one
	int otherRigidNum;
	MatVelocityField *rmvf=cvfm->GetRigidMaterialField(&otherRigidNum);
	if(rmvf==NULL) return;
	if(otherRigidNum!=rigidFieldNum)
		throw MPMTermination("Two different rigid materials on the same node","CrackVelocityFieldMulti::MaterialContact");
	
	// add number of rigid points and total points this crack velocity field
	numberRigidPoints+=rmvf->numberPoints;
	numberPoints+=rmvf->numberPoints;
	
	// add unscaled volume to this crack velocity field
	unscaledVolume+=cvfm->UnscaledVolumeRigid();
	
	// sum momentum, displacement, and mass grad (velocity is same) into material velocity field
	mvf[rigidFieldNum]->numberPoints+=rmvf->numberPoints;
	AddVector(&mvf[rigidFieldNum]->pk,&rmvf->pk);
	AddVector(&mvf[rigidFieldNum]->disp,&rmvf->disp);
	AddVector(mvf[rigidFieldNum]->massGrad,rmvf->massGrad);
}

// Copy rigid material from another velocity field (cvfm) to this field (creating if needed)
// This is only called if COMBINE_RIGID_MATERIALS is defined
void CrackVelocityFieldMulti::CopyRigidFrom(CrackVelocityFieldMulti *cvfm,int rigidFieldNum)
{	
	/*
	// create only if already in this field
	if(mvf[rigidFieldNum]==NULL) return;
	if(mvf[rigidFieldNum]->numberPoints==0) return;
	
	// current valiues
	int initialRigidPoints=mvf[rigidFieldNum]->numberPoints;
	unscaledVolume-=UnscaledVolumeRigid();
	*/
	
	// create material field if needed
	int initialRigidPoints=0;
	double initialRigidVolume=0.;
	if(mvf[rigidFieldNum]==NULL)
	{	mvf[rigidFieldNum]=new MatVelocityField(TRUE);
		if(mvf[rigidFieldNum]==NULL) throw CommonException("Memory error allocating material velocity field.",
													"CrackVelocityFieldMulti::CopyRigidFrom");
		numberMaterials++;					// just added a material to this crack velocity field
	}
	else
	{	initialRigidPoints=mvf[rigidFieldNum]->numberPoints;
		initialRigidVolume=UnscaledVolumeRigid();
		if(initialRigidPoints==0) numberMaterials++;
	}
	
	// reference to source field
	MatVelocityField *rmvf=cvfm->mvf[rigidFieldNum];
	
	// add number rigid points this crack velocity field
	numberRigidPoints+=rmvf->numberPoints-initialRigidPoints;
	numberPoints+=rmvf->numberPoints-initialRigidPoints;
	
	// add unscaled volume to this crack velocity field
	unscaledVolume+=cvfm->UnscaledVolumeRigid()-initialRigidVolume;
	
	// copy momentum, displacement, and mass grad (velocity is same) into material velocity field
	mvf[rigidFieldNum]->numberPoints=rmvf->numberPoints;
	CopyVector(&mvf[rigidFieldNum]->pk,&rmvf->pk);
	CopyVector(&mvf[rigidFieldNum]->vk,&rmvf->vk);
	CopyVector(&mvf[rigidFieldNum]->disp,&rmvf->disp);
	CopyVector(mvf[rigidFieldNum]->massGrad,rmvf->massGrad);
}
	
#pragma mark TASK 3 METHODS

// Add to fint spread out over the materials to each has same extra accerations = f/M
// only called to add interface force on a crack
void CrackVelocityFieldMulti::AddFintSpreadTask3(Vector *f)
{	int i;
	
	// special case for only one material
	if(numberMaterials==1)
	{	for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			{	AddVector(&mvf[i]->fint,f);
				break;
			}
		}
	}
	
	// more than one material, add to nonrigid materials only
	else
	{	double totMass=GetTotalMass();
		for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
				AddScaledVector(&mvf[i]->fint,f,mvf[i]->mass/totMass);
		}
	}
}

// Add to fext spread out over the materials to each has same extra accerations = f/M
// Only called for crack traction forces
void CrackVelocityFieldMulti::AddFextSpreadTask3(Vector *f)
{	int i;
	
	// special case for only one material
	if(numberMaterials==1)
	{	for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			{	AddVector(&mvf[i]->fext,f);
				break;
			}
		}
	}
	
	// more than one material
	else
	{	double totMass=GetTotalMass();
		for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
				AddScaledVector(&mvf[i]->fext,f,mvf[i]->mass/totMass);
		}
	}
}

// Calculate total force at a node from current values
// 		now m*a in g mm/sec^2
void CrackVelocityFieldMulti::CalcFtotTask3(double extDamping)
{	int i;
    for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			mvf[i]->CalcFtotTask3(extDamping);
	}
}

#pragma mark TASK 4 METHODS

// update momenta for this MPM step
//  pk(i+1) = pk(i) + ftot * dt
void CrackVelocityFieldMulti::UpdateMomentaOnField(double timestep)
{	// update momenta
	int i;
    for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			AddScaledVector(&mvf[i]->pk,&mvf[i]->ftot,timestep);
    }
}

#pragma mark TASK 6 METHODS

// zero momentum and displacement at a node for new calculations
// but can do the calculation for rigid particles here
void CrackVelocityFieldMulti::RezeroNodeTask6(double deltaTime)
{	int i;
    for(i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveField(mvf[i]))
		{	if(mvf[i]->mass>0.)
			{	ZeroVector(&mvf[i]->pk);
				ZeroVector(&mvf[i]->disp);
			}
			else
			{	// for rigid particles, can project displacement using current velocity because
				// particle mass is its volume
				// dnew = Sum (Vp*fpi*(d + v dt)) = dold + Sum (Vp*fpi*v*dt) = dold + pk*dt
				AddScaledVector(&mvf[i]->disp,&mvf[i]->pk,deltaTime);
			}
		}
    }
}

#pragma mark MATERIAL CONTACT

/* Called in multimaterial mode to check contact at nodes with multiple materials

	Input parameters:
		mvf[]->mass,pk,massGrad,disp (if contact by displacements)

	Output changes
		mvf[]->pk, vk (if one particle), ftot (if postUpdate is TRUE)
*/
void CrackVelocityFieldMulti::MaterialContact(int nodenum,int vfld,bool postUpdate,double deltime)
{
	// exit if no contact
	if(numberMaterials<=1) return;
	
	// get center of mass results and look out for rigid materials
	int i,j,rigidMat=-1;
	Vector Pc,dispc;
	ZeroVector(&Pc);
	double Mc=0.;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	double matMass=mvf[i]->mass;
			if(matMass>0.)
			{	// real material has mass
				AddVector(&Pc,&mvf[i]->pk);
				Mc+=mvf[i]->mass;
			}
			else if(rigidMat>0)
			{	// rigid material, but not allowed if already had another rigid material
				throw MPMTermination("Two different rigid materials in contact on the same node",
												"CrackVelocityFieldMulti::MaterialContact");
			}
			else
			{	// first rigid material at this node
				rigidMat=i;
			}
		}
	}
	
	// if exactly one rigid material, then special case for contact laws
	if(rigidMat>=0)
	{	RigidMaterialContact(rigidMat,nodenum,vfld,postUpdate,deltime);
		return;
	}
	
	// from here on all materials in contact are non-rigid
	if(contact.displacementCheck)
	{	dispc=GetCMDisplacement();
		ScaleVector(&dispc,1./Mc);
	}
	
	// loop over each material
	for(i=0;i<maxMaterialFields;i++)
    {	if(!MatVelocityField::ActiveField(mvf[i])) continue;
		
		// some variables
		Vector norm,delPi;
		double rho=MaterialBase::GetMVFRho(i);				// in g/mm^3
		double dotn,massi=mvf[i]->mass,massRatio=massi/Mc;
		
		// First determine contact law from other material with most volume
		double maxOtherMaterialVolume=0.,rhoj,rhopaired=rho;
		int ipaired=0;
		for(j=0;j<maxMaterialFields;j++)
		{	if(j==i || !MatVelocityField::ActiveField(mvf[j])) continue;
			rhoj=MaterialBase::GetMVFRho(j);				// in g/mm^3
			double matUnscaledVolume=mvf[j]->mass/rhoj;
			if(matUnscaledVolume>maxOtherMaterialVolume)
			{	maxOtherMaterialVolume=matUnscaledVolume;
				ipaired=j;
				rhopaired=rhoj;
			}
		}
		// problem if ipaired not found, but it will be found
		int maxContactLaw=contact.GetMaterialContactLaw(i,ipaired);
		double maxFriction=contact.GetMaterialFriction(i,ipaired);
		
		if(maxContactLaw!=NOCONTACT)
		{	// check nodal volume
			if(unscaledVolume/mpmgrid.GetCellVolume()<contact.materialContactVmin) continue;
		
			// ignore very small mass nodes
			if(massRatio<1.e-6 || massRatio>0.999999) continue;
			
			// find -mi(vi-vc) = (ma/mc)pc-pi or momentum change to match ctr of mass momentum
			CopyScaleVector(&delPi,&mvf[i]->pk,-1.);
			AddScaledVector(&delPi,&Pc,massRatio);

			// Get normal vector by various options
			switch(contact.materialNormalMethod)
			{	case MAXIMUM_VOLUME_GRADIENT:
				{	// Use mat with largest magnitude volume gradient
					Vector normi,normj;
					nd[nodenum]->GetMassGradient(vfld,i,&normi,1.);
					nd[nodenum]->GetMassGradient(vfld,ipaired,&normj,-1.);
					
					// compare magnitude of volume gradients
					double magi=sqrt(DotVectors(&normi,&normi));
					double magj=sqrt(DotVectors(&normj,&normj));
					if(magi/rho >= magj/rhopaired)
						CopyScaleVector(&norm,&normi,1./magi);		// use material i
					else
						CopyScaleVector(&norm,&normj,1./magj);		// use material j
					break;
				}
				case MAXIMUM_VOLUME:
					// Use mat with most volume
					if(massi/rho >= maxOtherMaterialVolume)
						nd[nodenum]->GetMassGradient(vfld,i,&norm,1.);
					else
						nd[nodenum]->GetMassGradient(vfld,ipaired,&norm,-1.);
					CopyScaleVector(&norm,&norm,1./sqrt(DotVectors(&norm,&norm)));
					break;

				default:
					break;
			}
			
			// Development code to try alternative methods
			if(fmobj->dflag[0]==1)
			{	// Use each material's own volume gradient
				nd[nodenum]->GetMassGradient(vfld,i,&norm,1.);
				ScaleVector(&norm,1./sqrt(DotVectors(&norm,&norm)));
			}
			else if(fmobj->dflag[0]==2)
			{	// get an average volume gradient
				Vector normj;
				nd[nodenum]->GetMassGradient(vfld,i,&norm,1.);
				nd[nodenum]->GetMassGradient(vfld,ipaired,&normj,-1.);
				AddVector(&norm,&normj);
				ScaleVector(&norm,1./sqrt(DotVectors(&norm,&norm)));
			}
			
			// get approach direction momentum form delPi.n (actual (vc-vi) = delPi/mi)
			dotn=DotVectors(&delPi,&norm);
			
			// With this check, any movement apart will be taken as noncontact
			// Also, frictional contact assumes dotn<0
			if(dotn>=0.) continue;
			
			// displacement check
			if(contact.displacementCheck)
			{	// get other mass and ignore if very small mass in other materials
				double scaleDisp=(Mc-massi)/Mc;
				if(scaleDisp<1.e-6) continue;
				
				// scale displacements to get delta = (Mc/(Mc-mi))*disp - (Mc/(Mc-mi))*(mvf[i]->disp/mi)
				Vector dispcScaled,dispi;
				CopyScaleVector(&dispcScaled,&dispc,1./scaleDisp);
				scaleDisp*=massi;
				CopyScaleVector(&dispi,&mvf[i]->disp,1./scaleDisp);
				
				// to get normal velocity delta v = (Mc/(Mc-mi)) (delta p/mi)
				double dvel = dotn/scaleDisp;
				
				// check for contact
				if(contact.MaterialContact(&dispi,&dispcScaled,&norm,dvel,postUpdate,deltime)==SEPARATED) continue;
			}
		}
		else
		{	// for no contact rule only get single velocity field conditions
			CopyScaleVector(&delPi,&mvf[i]->pk,-1.);
			AddScaledVector(&delPi,&Pc,massRatio);
		}
		
		// the material is in contact
		Vector tang;
		double dott;
		
		switch(maxContactLaw)
		{	case STICK:
			case NOCONTACT:
				break;
				
			case FRICTIONLESS:
				CopyScaleVector(&delPi,&norm,dotn);
				break;
				
			case FRICTIONAL:
				CopyVector(&tang,&delPi);
				AddScaledVector(&tang,&norm,-dotn);
				dott=sqrt(DotVectors(&tang,&tang));
				if(!DbleEqual(dott,0.))
				{	ScaleVector(&tang,1./dott);
					dott=DotVectors(&delPi,&tang);
					if(dott<0.)
					{	ScaleVector(&tang,-1.);
						dott=-dott;
					}
					if(dott>-maxFriction*dotn)
					{	AddScaledVector(&norm,&tang,-maxFriction);
						CopyScaleVector(&delPi,&norm,dotn);
					}
				}
				break;
				
			default:
				break;
		}
		
		// on post update contact, do not change nodes with boundary conditions
		unsigned char fixedDirection=nd[nodenum]->fixedDirection;
		if(postUpdate && fixedDirection)
		{	if(fixedDirection&X_DIRECTION) delPi.x=0.;
			if(fixedDirection&Y_DIRECTION) delPi.y=0.;
			if(fixedDirection&Z_DIRECTION) delPi.z=0.;
		}
		
		// change momenta
		mvf[i]->ChangeMatMomentum(&delPi,postUpdate,deltime);
		
		// special case two materials for efficiency (and if both will find normal the same way)
		if(numberMaterials==2 && fmobj->dflag[0]!=1)
		{	mvf[ipaired]->ChangeMatMomentum(ScaleVector(&delPi,-1.),postUpdate,deltime);
			break;
		}
	}
}

// Called in multimaterial mode to check contact at nodes with multiple materials and here
// means exactly one is a rigid material
//	(no rigid materials handled in MaterialContact(), two rigid materials is an error)
void CrackVelocityFieldMulti::RigidMaterialContact(int rigidFld,int nodenum,int vfld,bool postUpdate,double deltime)
{
	// get rigid material volume by subtracting other materials from the total unscaled volume
	double rho,rigidVolume=UnscaledVolumeRigid();
	int i;
	
	// loop over each material (skipping the one rigid material
	for(i=0;i<maxMaterialFields;i++)
    {	if(!MatVelocityField::ActiveField(mvf[i]) || i==rigidFld) continue;
		
		// First determine contact law with rigid material
		int maxContactLaw=contact.GetMaterialContactLaw(i,rigidFld);
		double maxFriction=contact.GetMaterialFriction(i,rigidFld);
		
		// NOCONTACT with rigid materials means to ignore contact, which it not revert ot single velocity field
		if(maxContactLaw==NOCONTACT) continue;
		
		// some variables
		Vector norm,delPi;
		rho=MaterialBase::GetMVFRho(i);				// in g/mm^3
		double dotn,massi=mvf[i]->mass;
		
		// check nodal volume
		if(unscaledVolume/mpmgrid.GetCellVolume()<contact.materialContactVmin) continue;
		
		// ignore very small interactions
		double volRatio=massi/rho/unscaledVolume;
		if(volRatio<1.e-6 || volRatio>0.999999) continue;
		
		// find -mi(vi-vr) = mi*vr-pi, which is change in momentum to match the rigid particle velocity
		CopyScaleVector(&delPi,&mvf[i]->pk,-1.);
		AddScaledVector(&delPi,&mvf[rigidFld]->vk,massi);
		
		// Get normal vector by various options
		switch(contact.materialNormalMethod)
		{	case MAXIMUM_VOLUME_GRADIENT:
			{	// Use mat with largest magnitude volume gradient
				Vector normi,normj;
				nd[nodenum]->GetMassGradient(vfld,i,&normi,1.);
				nd[nodenum]->GetMassGradient(vfld,rigidFld,&normj,-1.);
				
				// compare square of volume gradients (the bias has been squared)
				double magi=DotVectors(&normi,&normi);
				double magj=DotVectors(&normj,&normj);			// already a volume gradient
				if(magi/(rho*rho) >= contact.rigidGradientBias*magj)
					CopyScaleVector(&norm,&normi,1./sqrt(magi));		// use non-rigid material
				else
					CopyScaleVector(&norm,&normj,1./sqrt(magj));		// use rigid material
				break;
			}
			case MAXIMUM_VOLUME:
				// Use mat with most volume
				if(massi/rho >= rigidVolume)
					nd[nodenum]->GetMassGradient(vfld,i,&norm,1.);
				else
					nd[nodenum]->GetMassGradient(vfld,rigidFld,&norm,-1.);
				CopyScaleVector(&norm,&norm,1./sqrt(DotVectors(&norm,&norm)));
				break;
			/*
			case EACH_MATERIALS_MASS_GRADIENT:
				// Use each mat as is
				nd[nodenum]->GetMassGradient(vfld,i,&norm,1.);
				CopyScaleVector(&norm,&norm,1./sqrt(DotVectors(&norm,&norm)));
				break;
			 
			case AVERAGE_MAT_VOLUME_GRADIENTS:
				// Take an average of the two volume gradients
				break;
			*/
			default:
				break;
		}			
		
		// Development code to try alternative methods to get rigid contact normal
		if(fmobj->dflag[0]==1)
		{	// Use each material's own volume gradient
			nd[nodenum]->GetMassGradient(vfld,i,&norm,1.);
			ScaleVector(&norm,1./sqrt(DotVectors(&norm,&norm)));
		}
		else if(fmobj->dflag[0]==2)
		{	// get an average volume gradients (not same as averaging mass gradients with non-rigid materials)
			Vector normj;
			nd[nodenum]->GetMassGradient(vfld,i,&norm,1./rho);
			nd[nodenum]->GetMassGradient(vfld,rigidFld,&normj,-1.);
			AddVector(&norm,&normj);
			ScaleVector(&norm,1./sqrt(DotVectors(&norm,&norm)));
		}
		else if(fmobj->dflag[0]==3)
		{	// try special normals for cutting simulation
			Vector dispi;
			CopyScaleVector(&dispi,&mvf[i]->disp,1./massi);
			if(dispi.y<=1.90)
			{	norm.x=0.;
				norm.y=1.;
			}
			else
			{	norm.x=1./sqrt(2.);
				norm.y=-1./sqrt(2.);
			}
			norm.z=0.;
		}
		else if(fmobj->dflag[0]==4)
		{	// try special normals for cutting simulation with rake angle into dflag[1]
			Vector dispi;
			CopyScaleVector(&dispi,&mvf[i]->disp,1./massi);
			if(dispi.y<=0.0)
			{	norm.x=0.;
				norm.y=1.;
			}
			else
			{	double radAngle=(double)fmobj->dflag[1]*PI_CONSTANT/180.;
				norm.x=cos(radAngle);
				norm.y=-sin(radAngle);
			}
			norm.z=0.;
		}
		
		// get approach direction momentum form delPi.n (actual (vr-vi).n = delPi.n/mi)
		dotn=DotVectors(&delPi,&norm);
		
		// With this check, any movement apart will be taken as noncontact
		// Also, frictional contact assumes dotn<0
		if(dotn>=0.) continue;
		
		// displacement check
		if(contact.displacementCheck)
		{	// rigid material displacement was scaled by volume, while non-rigid was weighted by mass
			Vector rigidDisp,dispi;
			CopyScaleVector(&rigidDisp,&mvf[rigidFld]->disp,1./rigidVolume);
			CopyScaleVector(&dispi,&mvf[i]->disp,1./massi);
			
			// convert dotn to velocity of approach
			double dvel = dotn/massi;
			
			// check for contact
			if(contact.MaterialContact(&dispi,&rigidDisp,&norm,dvel,postUpdate,deltime)==SEPARATED) continue;
		}
		
		// the material is in contact
		Vector tang;
		double dott;
		
		switch(maxContactLaw)
		{	case STICK:
			case NOCONTACT:
				break;
				
			case FRICTIONLESS:
				CopyScaleVector(&delPi,&norm,dotn);
				break;
				
			case FRICTIONAL:
				CopyVector(&tang,&delPi);
				AddScaledVector(&tang,&norm,-dotn);
				dott=sqrt(DotVectors(&tang,&tang));
				if(!DbleEqual(dott,0.))
				{	ScaleVector(&tang,1./dott);
					dott=DotVectors(&delPi,&tang);
					if(dott<0.)
					{	ScaleVector(&tang,-1.);
						dott=-dott;
					}
					if(dott>-maxFriction*dotn)
					{	AddScaledVector(&norm,&tang,-maxFriction);
						CopyScaleVector(&delPi,&norm,dotn);
					}
				}
				break;
				
			default:
				break;
		}
		
		// on post update contact, do not change nodes with boundary conditions
		unsigned char fixedDirection=nd[nodenum]->fixedDirection;
		if(postUpdate && fixedDirection)
		{	if(fixedDirection&X_DIRECTION) delPi.x=0.;
			if(fixedDirection&Y_DIRECTION) delPi.y=0.;
			if(fixedDirection&Z_DIRECTION) delPi.z=0.;
		}
		
		// change momenta
		mvf[i]->ChangeMatMomentum(&delPi,postUpdate,deltime);
		
		// store contact force
		mvf[rigidFld]->AddContactForce(&delPi);
	}
}

	// retrieve mass gradient
void CrackVelocityFieldMulti::GetMassGradient(int matfld,Vector *grad,double scale) { CopyScaleVector(grad,mvf[matfld]->massGrad,scale); }

#pragma mark VELOCITY METHODS

// Calculate velocity at a node from current momentum and mass matrix in all velocity fields
void CrackVelocityFieldMulti::CalcVelocityForStrainUpdate(void)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			mvf[i]->CalcVelocityForStrainUpdate();
	}
}

#pragma mark BOUNDARY CONDITIONS

// zero x moment and velocity
void CrackVelocityFieldMulti::SetXMomVel(void)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	mvf[i]->pk.x=0.;
			mvf[i]->vk.x=0.;
		}
	}
}

// zero y moment and velocity
void CrackVelocityFieldMulti::SetYMomVel(void)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	mvf[i]->pk.y=0.;
			mvf[i]->vk.y=0.;
		}
	}
}

// zero z moment and velocity
void CrackVelocityFieldMulti::SetZMomVel(void)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	mvf[i]->pk.z=0.;
			mvf[i]->vk.z=0.;
		}
	}
}

// zero momentum in direction (cos(angle), -sin(angle)) or vector rotated from postive x axis
// by clockwise angle. The desired vector is (p.t)t where t is unit vector normal to
// skew direction and here t = (sin(angle), cos(angle))
void CrackVelocityFieldMulti::SetSkewMomVel(double angle)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	double c=cos(angle),s=sin(angle);
			Vector *npk=&mvf[i]->pk;
			double momx=npk->x*s*s + npk->y*c*s;
			double momy=npk->x*c*s + npk->y*c*c;
			npk->x=momx;
			npk->y=momy;
			mvf[i]->vk.x=momx/mvf[i]->mass;
			mvf[i]->vk.y=momy/mvf[i]->mass;
		}
	}
}

// add to x momentum and velocity from BCs
void CrackVelocityFieldMulti::AddXMomVel(double vx)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	mvf[i]->pk.x+=mvf[i]->mass*vx;
			mvf[i]->vk.x+=vx;
		}
	}
}

// add to y momentum and velocity from BCs
void CrackVelocityFieldMulti::AddYMomVel(double vy)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	mvf[i]->pk.y+=mvf[i]->mass*vy;
			mvf[i]->vk.y+=vy;
		}
	}
}

// add to z momentum and velocity from BCs
void CrackVelocityFieldMulti::AddZMomVel(double vz)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	mvf[i]->pk.z+=mvf[i]->mass*vz;
			mvf[i]->vk.z+=vz;
		}
	}
}

// Add velocity in direction (cos(angle), -sin(angle)) or vector rotated from postive x axis
// by clockwise angle.
void CrackVelocityFieldMulti::AddSkewMomVel(double vel,double angle)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	double velx=cos(angle)*vel;
			double vely=-sin(angle)*vel;
			mvf[i]->pk.x+=mvf[i]->mass*velx;
			mvf[i]->vk.x+=velx;
			mvf[i]->pk.y+=mvf[i]->mass*vely;
			mvf[i]->vk.y+=vely;
		}
	}
}

// set x force to -p(interpolated)/time such that updated momentum
//    of pk.x + deltime*ftot.x will be zero
void CrackVelocityFieldMulti::SetXFtot(double deltime)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			mvf[i]->ftot.x=-mvf[i]->pk.x/deltime;
	}
}

// set y force to -p(interpolated)/time such that updated momentum
//    of pk.y + deltime*ftot.y will be zero
void CrackVelocityFieldMulti::SetYFtot(double deltime)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			mvf[i]->ftot.y=-mvf[i]->pk.y/deltime;
	}
}

// set z force to -p(interpolated)/time such that updated momentum
//    of pk.z + deltime*ftot.z will be zero
void CrackVelocityFieldMulti::SetZFtot(double deltime)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			mvf[i]->ftot.z=-mvf[i]->pk.z/deltime;
	}
}

// Change current ftot such that updated momentum of (pk.x + deltime*ftot.x, pk.y + deltime*ftot.y) will have zero
// momentum in the (cos(theta), -sin(angle)) direction (or direction clockwise from positive x axis by angle).
// Superpose force to induce zero momentum in the skew direction:
//			f dt = (p.t)t - p(interpolated)
// where t = (sin(angle),cos(angle)) is tangential vector, with the existing component of total force
// in the t direction or (f.t)t
void CrackVelocityFieldMulti::SetSkewFtot(double deltime,double angle)
{	
	 int i;
	 for(i=0;i<maxCrackFields;i++)
	 {   if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	double c=cos(angle),s=sin(angle);
			double dfxdt=-mvf[i]->pk.x*c*c + mvf[i]->pk.y*c*s;	// to get zero skew momentum
			double dfydt=mvf[i]->pk.x*c*s - mvf[i]->pk.y*s*s;
			double fx=mvf[i]->ftot.x*s*s + mvf[i]->ftot.y*c*s;	// f normal to skew direction
			double fy=mvf[i]->ftot.x*c*s + mvf[i]->ftot.y*c*c;
			mvf[i]->ftot.x=fx + dfxdt/deltime;
			mvf[i]->ftot.y=fy + dfydt/deltime;
		}
	 }
}

// add to x force such that updated momentum will be mass*velocity
void CrackVelocityFieldMulti::AddXFtot(double deltime,double velx)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			mvf[i]->ftot.x+=mvf[i]->mass*velx/deltime;
	}
}

// add to y force such that updated momentum will be mass*velocity
void CrackVelocityFieldMulti::AddYFtot(double deltime,double vely)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			mvf[i]->ftot.y+=mvf[i]->mass*vely/deltime;
	}
}

// add to z force such that updated momentum will be mass*velocity
void CrackVelocityFieldMulti::AddZFtot(double deltime,double velz)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			mvf[i]->ftot.z+=mvf[i]->mass*velz/deltime;
	}
}

// set force in the skew direction (cos(angle),-sin(angle)), or direction rotated from postive
// x axis by clockwise angle.
void CrackVelocityFieldMulti::AddSkewFtot(double deltime,double vel,double angle)
{	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	double velx=cos(angle)*vel;
			double vely=-sin(angle)*vel;
			mvf[i]->ftot.x+=mvf[i]->mass*velx/deltime;
			mvf[i]->ftot.y+=mvf[i]->mass*vely/deltime;
		}
	}
}

#pragma mark ACCESSORS

// total number of non-rigid points
int CrackVelocityFieldMulti::GetNumberPointsNonrigid(void) { return numberPoints-numberRigidPoints; }

// location for crack in this field
// total mass all velocity fields (rigid particles have zero mass)
double CrackVelocityFieldMulti::GetTotalMass(void)
{	int i;
	double mass=0;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			mass+=mvf[i]->mass;
	}
	return mass;
}

// get one mass (rigid particles have zero mass
double CrackVelocityFieldMulti::GetMass(int matfld)
{	if(MatVelocityField::ActiveNonrigidField(mvf[matfld]))
		return mvf[matfld]->mass;
	else
		return 0.;
}

// get center of mass momentum for all nonrigid material fields in this crack velocity field
Vector CrackVelocityFieldMulti::GetCMatMomentum(void)
{	Vector pk;
	ZeroVector(&pk);
	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			AddVector(&pk,&mvf[i]->pk);
	}
	return pk;
}

// get center of mass displacement (actually sum of displacement*mass so displacement is vector/total mass)
// Includes only non-rigid materials
Vector CrackVelocityFieldMulti::GetCMDisplacement(void)
{	Vector dk;
	ZeroVector(&dk);
	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			AddVector(&dk,&mvf[i]->disp);
	}
	return dk;
}

// get center of mass force for all material fields in this crack velocity field
Vector CrackVelocityFieldMulti::GetCMatFtot(void)
{	Vector fk;
	ZeroVector(&fk);
	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			AddVector(&fk,&mvf[i]->ftot);
	}
	return fk;
}

// add contact force on rigid material to the input vector
void CrackVelocityFieldMulti::SumAndClearRigidContactForces(Vector *fcontact,bool clearForces)
{	int rigidFieldNum;
	MatVelocityField *rigidField=GetRigidMaterialField(&rigidFieldNum);
	if(rigidField!=NULL)
	{	AddVector(fcontact,&rigidField->ftot);
		if(clearForces) ZeroVector(&rigidField->ftot);
	}
}

// get first active rigid field or return NULL. Also return number in rigidFieldNum
// This is only called if COMBINE_RIGID_MATERIALS is defined
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

// get rigid material volume by subtracting other materials from the total unscaled volume
double CrackVelocityFieldMulti::UnscaledVolumeNonrigid(void)
{	// total volume if no rigid particles
	if(numberRigidPoints<=0) return unscaledVolume;
	
	// sum each nonrigid material
	double rho,nonrigidVolume=0.0;
	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	rho=MaterialBase::GetMVFRho(i);				// in g/mm^3
			nonrigidVolume+=(mvf[i]->mass/rho);
		}
	}
	return nonrigidVolume;
}

// get rigid material volume by subtracting nonrigid materials from the total unscaled volume
// assumes at most one rigid field
double CrackVelocityFieldMulti::UnscaledVolumeRigid(void)
{	// total volume if all rigid particles
	if(numberRigidPoints==numberPoints) return unscaledVolume;
	return unscaledVolume-UnscaledVolumeNonrigid();
}

/* in response to crack contact, change moment by changing velocity of all 
	nonrigid materials the same amount
 
   Change velocity by dP/M, where M is total mass
   Material i velocity becomes vi = pi/mi + dP/M
   Material i momentum change is mi vi = pi + mi dP/M
*/
void CrackVelocityFieldMulti::ChangeMomentum(Vector *delP,bool postUpdate,double deltime)
{
	int i;
	
	// special case for only one material (and it must be nonrigid
	if(numberMaterials==1)
	{	for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveField(mvf[i]))
			{	mvf[i]->ChangeMatMomentum(delP,postUpdate,deltime);
				break;
			}
		}
	}
	
	// more than one material
	else
	{	Vector partialDelP;
		double totMass=GetTotalMass();
		for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
				mvf[i]->ChangeMatMomentum(CopyScaleVector(&partialDelP, delP, mvf[i]->mass/totMass),postUpdate,deltime);
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
void CrackVelocityFieldMulti::Describe(void)
{	CrackVelocityField::Describe();
	cout << "#     multimaterial: nmat= " << numberMaterials << " nrigidpts= " << numberRigidPoints << endl;
	int i;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			mvf[i]->Describe();
	}
}


	