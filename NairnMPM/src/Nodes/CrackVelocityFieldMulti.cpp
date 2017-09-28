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
		{
			delete mvf[i];
		}
	}
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

// match material velocity fields on ghost node to those on a real node
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

// Copy rigid material from field [0] to this field (creating if needed)
// (note: can search "combining_rigid" to see all places in code affected by combining rigid particles)
// throws CommonException on material velocity allocation memory error
void CrackVelocityFieldMulti::CopyRigidFrom(MatVelocityField *rmvf,int rigidFieldNum)
{	
	// create material field in this crack velocity field if needed, otherwise, just be sure it is zeroed.
	if(mvf[rigidFieldNum]==NULL)
	{	mvf[rigidFieldNum] = new MatVelocityField(RIGID_FIELD_BIT);
		if(mvf[rigidFieldNum]==NULL) throw CommonException("Memory error allocating material velocity field.",
													"CrackVelocityFieldMulti::CopyRigidFrom");
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
	mvf[rigidFieldNum]->SetVelocity(&rvel);
    
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
// Only addes force to fields that see cracks
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
	else
	{	double totMass=GetTotalMass(true);
		for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
				mvf[i]->AddFtotScaled(f,mvf[i]->mass/totMass);
		}
	}
}

// Add gravity and body force at a node in g mm/sec^2 (non rigid only)
// For materials that ignore cracks, only add to source MVF in field [0]
void CrackVelocityFieldMulti::AddGravityAndBodyForceTask3(Vector *gridBodyForce)
{	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
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
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			mvf[i]->UpdateMomentum(timestep);
    }
}

#pragma mark MATERIAL CONTACT

/* Called in multimaterial mode to check contact at nodes with multiple materials

	Input parameters:
		mvf[]->mass,pk,volumeGrad,disp (if contact by displacements)

	VEL1 no longer needs to change vk in single particle node
	Output changes are only allowed on this node (to be thread safe for parallel)
		changes on mvf[]: pk, vk (if one particle), ftot (if postUpdate is TRUE)
 
	ndptr is parent node to this crack velocity field
 
	On first call in time step, first and last are pointers to Cracknode * because those
		objects are created for later interface calculations
	postUpdate is TRUE when called between momentum update and particle update and otherwise is false
	throws std::bad_alloc
*/
void CrackVelocityFieldMulti::MaterialContactOnCVF(NodalPoint *ndptr,double deltime,int callType)
{
	// exit if no contact
	if(numberMaterials<=1) return;
	
	// get center of mass results and look out for rigid materials
	int i,j,rigidMat=-1;
    bool multiRigid = false;
	Vector Pc,dispc;
	ZeroVector(&Pc);
	double Mc=0.;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	if(!mvf[i]->IsRigidField())
			{	// real material
				AddVector(&Pc,&mvf[i]->pk);
				Mc += mvf[i]->mass;
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
    {   //if(callType!=UPDATE_MOMENTUM_CALL) return;
		RigidMaterialContactOnCVF(rigidMat,multiRigid,ndptr,deltime,callType);
		return;
	}
	
	// from here on all materials in contact are non-rigid
	AdjustForSymmetry(ndptr,&Pc,false);
	
	// will probably need center of mass if doing displacement check, if imperfect interface,
	// or if transport tasks are doing contact. Because this is amost always, it is always
	// found
	dispc = GetCMDisplacement(ndptr,false);
	ScaleVector(&dispc,1./Mc);
	AdjustForSymmetry(ndptr,&dispc,false);
	
	// total volume for contact
	double totalVolume = GetVolumeTotal(ndptr);
    
	// loop over each material
	bool hasInterfaceEnergy = false;
    double qrate = 0.;
	bool postUpdate = callType != MASS_MOMENTUM_CALL;
	for(i=0;i<maxMaterialFields;i++)
	{	// continue if not an active field
		if(!MatVelocityField::ActiveField(mvf[i])) continue;
		
		// some variables
		Vector norm,delPi,delta,tangDel;
		Vector dispa,dispb;				// extrapolated position/mass for this material (a==i) and other material b
		double deln=0.;					// delta.n corrected for extrapolation method
		double delt=0.;					// delta.t
		double dotn,massi=mvf[i]->mass,massRatio=massi/Mc;
		double mred = (Mc-massi)/Mc;	// final mred = massi*(Mc-massi)/Mc, but found later when needed
        double voli=mvf[i]->GetContactVolume();
		double contactArea = -1.,mredDE = -1.,contactGridN = 1.;
		
		// Determine contact law from other material with most volume
        // Find lumped volume and volume gradient of all other non-rigid materials
		Vector otherGrad,normj;
		ZeroVector(&otherGrad);
		double otherVolume = totalVolume - mvf[i]->GetContactVolume();
		double maxOtherMaterialVolume = 0.;
		int ipaired = -1;
		for(j=0;j<maxMaterialFields;j++)
		{	if(j==i || !MatVelocityField::ActiveField(mvf[j])) continue;
			double matVolume = mvf[j]->GetContactVolume();
			if(matVolume>maxOtherMaterialVolume)
			{	maxOtherMaterialVolume = matVolume;
				ipaired = j;
			}
			// Finding  - Sum grad V_j
			GetVolumeGradient(j,ndptr,&normj,-1.);
			AddVector(&otherGrad,&normj);
		}
		
		// problem if ipaired not found, but it will always be found
		if(ipaired<0) continue;
		ContactLaw *theContactLaw = contact.GetMaterialContactLaw(i,ipaired);
		
		// Determine if in contact
        bool inContact = false;
		bool comContact = false;
        
        // If needed, determine if the surfaces are in contact
        if(theContactLaw->IgnoreContact() || comContact)
        {   // find -mi(vi-vc) = (ma/mc)pc-pi or momentum change to match ctr of mass momentum
            CopyScaleVector(&delPi,&mvf[i]->pk,-1.);
			AdjustForSymmetry(ndptr, &delPi, false);
            AddScaledVector(&delPi,&Pc,massRatio);
            inContact = true;
			comContact = true;
        }
        
        else
		{	//----------------------------------------------------------------------------------
			// 1. check nodal volume (this is turned off by setting the materialContactVmin to zero)
			//    (warning: 2D must set grid thickness if it is not 1)
			if(contact.materialContactVmin>0.)
			{	if(totalVolume/mpmgrid.GetCellVolume(ndptr)<contact.materialContactVmin) break;
			}
		
			//----------------------------------------------------------------------------------
			// 2. ignore very small mass nodes - may not be needed
			if(massRatio<1.e-6 || massRatio>0.999999) continue;
		
			//----------------------------------------------------------------------------------
			// 3. go through contact conditions; break if not in contact or
			//    set inContact to true and break if is in contact. Note that
			//    imperfect interfaces will always proceed to calculations, but it
			//    needs normal vector, displacements, and needs to know if in contact (to know
			//    whether to use Dnc or Dnt for normal stiffness). This never break without
			//	  these things calculated
			while(true)
			{	// find -mi(vi-vc) = (ma/mc)pc-pi or momentum change to match ctr of mass momentum
                CopyScaleVector(&delPi,&mvf[i]->pk,-1.);
				AdjustForSymmetry(ndptr, &delPi, false);
                AddScaledVector(&delPi,&Pc,massRatio);

				//----------------------------------------------------------------------------------
				// 3A. Get normal vector by various options
                switch(contact.materialNormalMethod)
                {	case MAXIMUM_VOLUME_GRADIENT:
                    {	// Use mat with largest magnitude volume gradient
                        Vector normi;
                        GetVolumeGradient(i,ndptr,&normi,1.);
                        
                        // compare magnitude of volume gradients
                        double magi=sqrt(DotVectors(&normi,&normi));
                        double magj=sqrt(DotVectors(&otherGrad,&otherGrad));
                        if(magi >= magj)
                            CopyScaleVector(&norm,&normi,1./magi);			// use material i
                        else
                            CopyScaleVector(&norm,&otherGrad,1./magj);		// use material j
                        break;
                    }
                        
                    case MAXIMUM_VOLUME:
						// Use mat with most volume
                        if(voli >= otherVolume)
                            GetVolumeGradient(i,ndptr,&norm,1.);
                        else
                            CopyVector(&norm,&otherGrad);
                        CopyScaleVector(&norm,&norm,1./sqrt(DotVectors(&norm,&norm)));
						break;
                        
                    case AVERAGE_MAT_VOLUME_GRADIENTS:
                    {	// get mass gradients material i and for other material(s)
                        Vector normi;
                        GetVolumeGradient(i,ndptr,&normi,1.);
                        
                        // volume weighted mean of volume gradients
                        //  = (voli * grad voli + otherVolume * grad otherVolume)/(total volume)
                        CopyScaleVector(&norm,&normi,voli);
                        AddScaledVector(&norm,&otherGrad,otherVolume);
                        
                        // normalize
                        double magi=sqrt(DotVectors(&norm,&norm));
                        ScaleVector(&norm,1./magi);
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
				
				//----------------------------------------------------------------------------------
				// 3A+. Various developer flag hacks in OSParticulas here

				//----------------------------------------------------------------------------------
				// 3B. Find relative velocity change. <0 is separating and always assumed then
				// to not be in contact.
				
				// get approach direction momentum from delPi.n (actual (vc-vi) = delPi/mi)
				dotn=DotVectors(&delPi,&norm);
				
				// Before exiting on dotn>=0, see if displacement calculations are needed
				//		(or for now, just do them)
				// We must get displacement if:
				//   a. A displacement check is needed
				//   b. contact heating is enabled
				//   c. contact law needs area and it will continue even if not in contact
				
				//----------------------------------------------------------------------------------
				// 3C. Do displacement calculations
				// Get COD and displacement of each material and final reduced mass
				// Let dispc = xiCOM, dispa = xia = mvf[i]->disp/massi, and mred = (Mc-Mi)/Mc
				// From my 2013 paper, COD = xib - xia = (Mc/(Mc-mi))(xiCOM - xia) = (dispc-xia)/mred
				
				// get dispa for this material
				CopyScaleVector(&dispa,&mvf[i]->disp,1./massi);
				AdjustForSymmetry(ndptr,&dispa,false);

				// get delta (dispc already adjusted for symmetry)
				CopyScaleVector(&delta,&dispc,1./mred);
				AddScaledVector(&delta,&dispa,-1./mred);
				
				// Final reduced mass = (mi*(Mc-mi))/Mc while until here was (Mc-mi)/Mc
				mred *= massi;
				
				// get dispb = delta + xia = xib-xia+xia = xib (already adjusted for symmetry)
				// for "other" material
				dispb = delta;
				AddVector(&dispb,&dispa);
				
				// get deln = delta.n corrected if simulation is extrapolating positiong
				deln = contact.MaterialSeparation(&delta,&dispa,&dispb,&norm,ndptr);
				
				//----------------------------------------------------------------------------------
				// 3B+. Now break not in contact if moving apart
				if(dotn>=0.) break;
				
				// ----------------------------------------------------------------------------
				// 3D. Displacement check:
                if(contact.displacementCheck)
				{	// on post update, adjust by normal velocity difference
					// to get normal velocity delta v = (Mc/(Mc-mi)) (delta p/mi)
					if(postUpdate)
					{	double dvel = dotn/mred;
						if(deln+dvel*deltime>=0.) break;
					}
					else
					{	// if current displacement positive then no contact
						if(deln>=0.) break;
					}
				}
				
				// ----------------------------------------------------------------------------
                // 3E. passed all tests, so in contact and break
                inContact = true;
                break;
            }
			
			// nan here means volume gradients zero with most common cause being gradient
			// normal to a symmetry plane
			if(norm.x!=norm.x || norm.y!=norm.y || norm.z!=norm.z)
			{	// no choice but to skip it
				continue;
			}
			
			// ----------------------------------------------------------------------------
			// 4. if contact is being used in transport tasks, OSParticuals only

			// ----------------------------------------------------------------------------
			// 5. if not in contact, not imperfect interface, and not adhesion (last two need to contine), then
			// done if only two materials (and the two use same normal) or continue to next material
			if(theContactLaw->ContactIsDone(inContact))
			{	if(numberMaterials==2 && contact.materialNormalMethod!=EACH_MATERIALS_MASS_GRADIENT) break;
				continue;
			}
		}
		
		// ----------------------------------------------------------------------------------
		// 6. Here inContact is usually true, but may be false for imperfet interface or for
		// anuy other contact law that want to continue (e.g., friction with adhesion)
		if(comContact)
		{	// ----------------------------------------------------------------------------------
			// 6A. use delPi found above and no further change needed, no other data needed either
		}
		else if(theContactLaw->IsFrictionalContact())
		{	// ----------------------------------------------------------------------------------
			// 6B. Handle frictional law
			
			// get area if needed this friction law
			if(contactArea<0. && theContactLaw->ContactLawNeedsContactArea())
			{	contactArea = GetContactArea(ndptr,voli,GetVolumeNonrigid(false)-voli,
											 &delta,&norm,&tangDel,&delt,&contactGridN);
			}
			
			// second order heating needs acceleration too
			// Find dF = (ma Fb - mb Fa)/Mc = (ma Fc - Mc Fa)/Mc
			// note that dacc = dF/mred
			bool getHeating = (callType==UPDATE_MOMENTUM_CALL) && ConductionTask::matContactHeating;
			Vector delFi;
			Vector *delFiPtr = NULL;
			if(getHeating)
			{	delFi = GetCMatFtot();
				ScaleVector(&delFi, massi);
				AddScaledVector(&delFi, mvf[i]->GetFtotPtr(), -Mc);
				ScaleVector(&delFi,1./Mc);
				delFiPtr = &delFi;
			}
			
			// get change in momentum, but false return means friction law want no change (i.e, continue as if not in contact)
			if(!theContactLaw->GetFrictionalDeltaMomentum(&delPi,&norm,dotn,&mredDE,mred,getHeating,
														  contactArea,inContact,deltime,delFiPtr))
			{	if(numberMaterials==2 && contact.materialNormalMethod!=EACH_MATERIALS_MASS_GRADIENT) break;
				continue;
			}
		}
		else
		{	// ----------------------------------------------------------------------------------
			// 6B. Handle imperfect interface
			
			// Interfaces always need contact area and tangDel and delt (find if not done already)
			if(contactArea<0.)
			{	contactArea = GetContactArea(ndptr,voli,GetVolumeNonrigid(false)-voli,
											 &delta,&norm,&tangDel,&delt,&contactGridN);
			}

			// get input force when needed and then get interface force and energy
			// Find delFi = (ma Fb - mb Fa)/Mc = (ma Fc - Mc Fa)/Mc = (ma Fc/Mc) - Fa (when needed)
			Vector fImp;
			if(callType==UPDATE_MOMENTUM_CALL)
			{	fImp = GetCMatFtot();
				ScaleVector(&fImp, massi/Mc);
				AddScaledVector(&fImp, mvf[i]->GetFtotPtr(), -1.);
			}
			else
				ZeroVector(&fImp);
			double rawEnergy;
			theContactLaw->GetInterfaceForces(&norm,&fImp,&rawEnergy,contactArea,&delPi,dotn,mred,
											  &tangDel,deln,delt,contactGridN);
			
#ifndef MANDMIMPINT
			if(callType==UPDATE_MOMENTUM_CALL)
			{	// add force (if any) to momentum change
				AddScaledVector(&delPi, &fImp, timestep);
				
				// Add interface energy. (Legacy units g-mm^2/sec^2 or multiply by 1e-9 to get J - kg-m^2/sec^2)
				NodalPoint::interfaceEnergy += rawEnergy;
				
				// but energy only once for this node
				hasInterfaceEnergy = true;
			}
#else
			if(callType==MASS_MOMENTUM_CALL)
			{	// add total force to material field (Legacy units microN)
				mvf[i]->AddFtot(&fImp);
				
				// decide if force balance can get other node too
				int iother = numberMaterials==2 && contact.materialNormalMethod!=EACH_MATERIALS_MASS_GRADIENT ? ipaired : -1 ;
				
				// add negative force to paired material (in a pair)
				if(iother>=0)
				{   ScaleVector(&fImp, -1.);
					mvf[iother]->AddFtot(&fImp);
				}
				
				// add interface energy if has paired material or if not added before
				if(iother>=0 || !hasInterfaceEnergy)
				{	// Add interface energy. (Legacy units g-mm^2/sec^2 or multiply by 1e-9 to get J - kg-m^2/sec^2)
					NodalPoint::interfaceEnergy += rawEnergy;
				}
				
				// has energy at least once
				hasInterfaceEnergy = true;
			}
#endif
		}
		
		// ----------------------------------------------------------------------------------
		// 7. Apply momentum change to the velocity fields
		
		// do not change nodes with boundary conditions
		ndptr->AdjustDelPiForBCs(&delPi);
		
		// change momenta
		mvf[i]->ChangeMatMomentum(&delPi,postUpdate,deltime);
		
		// special case two materials for efficiency (and if both will find normal the same way)
		if(numberMaterials==2 && contact.materialNormalMethod!=EACH_MATERIALS_MASS_GRADIENT)
		{	mvf[ipaired]->ChangeMatMomentum(ScaleVector(&delPi,-1.),postUpdate,deltime);
            
            // only true when in momentum update and conduction is on
            if(mredDE>0.)
			{	qrate += mredDE/mred;
            }

			break;
		}
        
        // friction when more than two materials (only true when in momentum update and conduction is on)
        if(mredDE>0.)
        {   qrate += mredDE/massi;
        }
	}
	
    // total friction, scale and absolute value
    if(callType==UPDATE_MOMENTUM_CALL && ConductionTask::matContactHeating && qrate!=0.)
    {	NodalPoint::frictionWork += qrate;
		conduction->AddFluxCondition(ndptr,fabs(qrate/deltime),true);
    }
}

// Called in multimaterial mode to check contact at nodes with multiple materials and here
// contact with one or more rigid matrerials (no rigid materials handled in MaterialContactOnCVF()
// throws std::bad_alloc
void CrackVelocityFieldMulti::RigidMaterialContactOnCVF(int rigidFld,bool multiRigid,NodalPoint *ndptr,double deltime,int callType)
{
    // multiRigid means more than one rigid material
	int i;
    double rigidVolume = mvf[rigidFld]->GetContactVolume();
    if(multiRigid)
    {   // this aborts the calculation
        //throw CommonException("Two different rigid materials in contact on the same node",
        //						"CrackVelocityFieldMulti::RigidMaterialContactOnCVF");
        
        // scheme 1, find node with the most contact volume
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
	
	// rigid material with position-dependent velocities may have mixed velocities
	Vector rvel;
	CopyScaleVector(&rvel, &mvf[rigidFld]->pk, 1./actualRigidVolume);
	// VEL1 special case for single particle node
	//if(mvf[rigidFld]->numberPoints==1)
    //    rvel = mvf[rigidFld]->GetVelocity();
	//else
	//	CopyScaleVector(&rvel, &mvf[rigidFld]->pk, 1./actualRigidVolume);
	AdjustForSymmetry(ndptr,&rvel,false);
	
	// loop over each material (skipping the one rigid material)
	bool postUpdate = callType != MASS_MOMENTUM_CALL;
	for(i=0;i<maxMaterialFields;i++)
    {	if(!MatVelocityField::ActiveField(mvf[i]) || mvf[i]->IsRigidField()) continue;
        
		// First determine contact law with rigid material
		ContactLaw *theContactLaw = contact.GetMaterialContactLaw(i,rigidFld);
        
        // get contact parameters
		Vector delta,tangDel;
		Vector dispRigid,dispi;				// displacement of rigid material and this material
		double deln=0.;						// delta.n corrected of extrapolation method
		double delt=0.;						// delta.t
        bool inContact = false;
 	
		// some variables
		Vector norm,delPi;
		double dotn,massi=mvf[i]->mass,voli=mvf[i]->GetContactVolume();
		double contactArea = -1.,mredDE = -1.,contactGridN = 1.;
        
        // find -mi(vi-vr) = - pi + mi*vr, which is change in momentum to match the rigid particle velocity
        CopyScaleVector(&delPi,&mvf[i]->pk,-1.);
		AdjustForSymmetry(ndptr,&delPi,false);
        AddScaledVector(&delPi,&rvel,massi);
		
        // if needed, determine if the surfaces are in contact
		bool comContact = false;
        if(theContactLaw->IgnoreContact())
        {   // will use -mi(vi-vr) = - pi + mi*vr to match rigid particle velocity
            inContact = true;
			comContact = true;
        }
        
        else
        {   //----------------------------------------------------------------------------------
            // 1. check nodal volume (this is turned off by setting the materialContactVmin to zero)
            //    (warning: 2D must set grid thickness if it is not 1)
			if(contact.materialContactVmin>0.)
            {	if(GetVolumeTotal(ndptr)/mpmgrid.GetCellVolume(ndptr)<contact.materialContactVmin) continue;
			}
            
			//----------------------------------------------------------------------------------
            // 2. ignore very small interactions
            double volRatio=voli/(rigidVolume+GetVolumeNonrigid(false));
            if(volRatio<.001 || volRatio>0.999) continue;
            //if(volRatio<1.e-6 || volRatio>0.999999) continue;
            
			//----------------------------------------------------------------------------------
            // 3. go through contact conditions; break if not in contact or
            //    set inContact to true and break if is in contact. Note that
            //    imperfect interfaces will always proceed to calculations, but it
            //    needs normal vector, interface displacements, and needs to know if in contact
			//	  (to know whether to use Dnc or Dnt for normal stiffness)
            while(true)
			{	//----------------------------------------------------------------------------------
				// 3A. Get normal vector by various options
                switch(contact.materialNormalMethod)
                {	case MAXIMUM_VOLUME_GRADIENT:
                    {	// Use mat with largest magnitude volume gradient
                        Vector normi,normj;
                        GetVolumeGradient(i,ndptr,&normi,1.);
                        GetVolumeGradient(rigidFld,ndptr,&normj,-1.);
                    
                        // compare square of volume gradients (the bias has been squared)
                        double magi=DotVectors(&normi,&normi);
                        double magj=DotVectors(&normj,&normj);
                        if(magi >= contact.rigidGradientBias*magj)				// squared here to avoid extra sqrt()
                            CopyScaleVector(&norm,&normi,1./sqrt(magi));		// use non-rigid material
                        else
                            CopyScaleVector(&norm,&normj,1./sqrt(magj));		// use rigid material
                        break;
                    }
                    
                    case MAXIMUM_VOLUME:
                        // Use mat with most volume (no rigid bias used)
                        if(voli >= rigidVolume)
                            GetVolumeGradient(i,ndptr,&norm,1.);
                        else
							GetVolumeGradient(rigidFld,ndptr,&norm,-1.);
                        CopyScaleVector(&norm,&norm,1./sqrt(DotVectors(&norm,&norm)));
						break;
                        
                    case AVERAGE_MAT_VOLUME_GRADIENTS:
                    {	// get volume-weighted mean of volume gradiants of material and rigid material
                        Vector normi,normj;
                        GetVolumeGradient(i,ndptr,&normi,1.);
                        GetVolumeGradient(rigidFld,ndptr,&normj,-1.);
                    
                        // volume weighted mean of volume gradients
                        //  = voli*normi + rigidVolume*normj (then normalized)
                        CopyScaleVector(&norm,&normi,voli);
                        AddScaledVector(&norm,&normj,rigidVolume);
                        double sumVolume = voli + rigidVolume;
                        double magi = DotVectors(&norm,&norm);
                        double magj = DotVectors(&normj,&normj);

                        // compare square of volume gradients (the bias has been squared)
                        if(magi/(sumVolume*sumVolume) >= contact.rigidGradientBias*magj)
                            ScaleVector(&norm,1./sqrt(magi));
                        else
                            CopyScaleVector(&norm,&normj,1./sqrt(magj));		// use rigid material
                        break;
                    }
                    
                    case EACH_MATERIALS_MASS_GRADIENT:
                        // Use non-rigid material's own gradient (no rigid bias used)
                        GetVolumeGradient(i,ndptr,&norm,1.);
                        CopyScaleVector(&norm,&norm,1./sqrt(DotVectors(&norm,&norm)));
                        break;

					case SPECIFIED_NORMAL:
						// use specified normal for all contact
						norm = contact.contactNormal;
						AdjustForSymmetry(ndptr,&norm,true);
						break;
                        
                    default:
                        break;
                }

				//----------------------------------------------------------------------------------
				// 3A+. Various developer flag hacks in OSParticulas are here
				
				//----------------------------------------------------------------------------------
				// 3B. Find relative velocity change. <0 is separating and always assume then
				// to not be in contact.
				
				// get approach direction momentum from delPi.n (actual (vc-vi) = delPi/mi)
				dotn=DotVectors(&delPi,&norm);
				
				// Before exiting on dotn>=0, see if displacement calculations are needed
				//		(or for now, just do them)
				// We must get displacement if:
				//   a. A displacement check is needed
				//   b. contact law needs area and it will continue even if not in contact
				
				//----------------------------------------------------------------------------------
				// 3C. Displacement calculations
				
				// rigid material displacement was scaled by volume, while non-rigid was weighted by mass
				CopyScaleVector(&dispRigid,&mvf[rigidFld]->disp,1./actualRigidVolume);
				AdjustForSymmetry(ndptr,&dispRigid,false);
				CopyScaleVector(&dispi,&mvf[i]->disp,1./massi);
				AdjustForSymmetry(ndptr,&dispi,false);
				delta = dispRigid;
				SubVector(&delta,&dispi);
				
				// get delta.n
				deln = contact.MaterialSeparation(&delta,&dispi,&dispRigid,&norm,ndptr);
				
				//----------------------------------------------------------------------------------
				// 3B+. Now break not in contact if moving apart
				if(dotn>=0.) break;
				
				//----------------------------------------------------------------------------------
                // 3D. Displacement check
                if(contact.displacementCheck)
                {	// on post update, adjust by normal velocity difference
					if(postUpdate)
					{	double dvel = dotn/massi;
						if(deln+dvel*deltime >= 0.) break;
					}
					else
					{	// if current displacement positive then no contact
						if(deln >= 0.) break;
					}
				}
                
                // passed all tests
                inContact = true;
                break;
            }
            
			// nan here means volume gradients zero with most common cause being gradient
			// normal to a symmetry plane
			if(norm.x!=norm.x || norm.y!=norm.y || norm.z!=norm.z)
			{	// no choice but to skip it
				continue;
			}
			
			// ----------------------------------------------------------------------------
			// 4. if not in contact, not imperfect interface, and not adhesion
			// (last two need to contine), then done now
			if(theContactLaw->ContactIsDone(inContact)) continue;
        }
		
		// ----------------------------------------------------------------------------------
		// 5. Here inContact is usually true, but may be false for imperfet interface or for
		// anuy other contact law that want to continue (e.g., friction with adhesion)
		if(comContact)
		{	// ----------------------------------------------------------------------------------
			// 5A. use delPi found above and no further change needed, no other data needed either
		}
		else if(theContactLaw->IsFrictionalContact())
		{	// ----------------------------------------------------------------------------------
			// 5B. Handle frictional law
			
			// may need area and contactGridN, but get contactGridN even if do not need area
			if(contactArea<0. && theContactLaw->ContactLawNeedsContactArea())
			{	contactArea = GetContactArea(ndptr,voli,rigidVolume,&delta,&norm,&tangDel,&delt,&contactGridN);
			}
			else
			{	Vector gridDist = mpmgrid.GetPerpendicularDistance(&norm,ndptr);
				contactGridN = gridDist.x;
			}
			
			// Second order heating needs to find delFi = -fi,a
			Vector delFi;
			Vector *delFiPtr = NULL;
			bool getHeating = (callType==UPDATE_MOMENTUM_CALL) && ConductionTask::matContactHeating;
			if(getHeating)
			{	delFi = mvf[i]->GetFtot();
				ScaleVector(&delFi, -1.);
				delFiPtr = &delFi;
			}
			
			// get change in momentum, but false return means friction law wants no change (i.e, continue as if not in contact)
			if(!theContactLaw->GetFrictionalDeltaMomentum(&delPi,&norm,dotn,&mredDE,massi,getHeating,
														  contactArea,inContact,deltime,delFiPtr))
				continue;
			
			// The following section is inspired Yang theses (U of W with Peter Arduino) which scaled
			// force (i.e. delPi) depending on distance from node to edge of the rigid material
			// Yang defined a rigid body. Here rigid body is collect of rigid particles
			double da;
			if(mpmgrid.GetContactByDisplacements())
			{	// need to convert displacement into distance ?
				da = -1.;
			}
			else
			{	// find extrapolated distance from rigid material to the node
				double pa = (dispRigid.x-ndptr->x)*norm.x + (dispRigid.y-ndptr->y)*norm.y + (dispRigid.z-ndptr->z)*norm.z;
				double r = mpmgrid.positionCutoff;
				if(r>0.)
					da = pa>0. ? pa/contactGridN - 0.5*r : pa/contactGridN + 0.5*r ;
				else
				{	r = -r;
					da = pa>0. ? 2.*pow(pa/(1.25*contactGridN),r) - 1. : 1 - 2.*pow(-pa/(1.25*contactGridN),r);
				}
			}
			
			// scale from 1 at 0 to 0 at 1
			double cutoff=0.5;
			if(da>cutoff)
			{	// rescale to zone from cutoff (0) to 1 (1)
				double oneMinusDa = 1. - fmin((da-cutoff)/(1.-cutoff),1.);
				
				//ScaleVector(&delPi,0.);							// truncate
				ScaleVector(&delPi,oneMinusDa);						// linear
				//ScaleVector(&delPi,oneMinusDa*oneMinusDa);		// square
				//ScaleVector(&delPi,sqrt(oneMinusDa));				// sqrt
				//ScaleVector(&delPi,pow(oneMinusDa,4));			// power
				
				// sigmoidal
				//double sc=12.,snorm=(1+exp(0.5*sc))/(1+exp(-0.5*sc));
				//double sigmoid = (1+exp(0.5*sc))/(1+exp(sc*(0.5-oneMinusDa)));
				//ScaleVector(&delPi,(sigmoid-1.)/(snorm-1.));
			}
		}
		else
		{	// ----------------------------------------------------------------------------------
			// 5B. Handle imperfect interface
			
			// Interfaces always need contact area (find if not done already)
			if(contactArea<0.)
			{	contactArea = GetContactArea(ndptr,voli,rigidVolume,&delta,&norm,&tangDel,&delt,&contactGridN);
			}

			// get input force when needed and then get interface force and energy
			// Find delFi = (ma Fb - mb Fa)/Mc = (ma Fc - Mc Fa)/Mc = (ma Fc/Mc) - Fa = -Fa (when needed)
			Vector fImp;
			ZeroVector(&fImp);
			if(callType==UPDATE_MOMENTUM_CALL)
				AddScaledVector(&fImp, mvf[i]->GetFtotPtr(), -1.);
			double rawEnergy;
			theContactLaw->GetInterfaceForces(&norm,&fImp,&rawEnergy,
													contactArea,&delPi,dotn,massi,&tangDel,deln,delt,contactGridN);
#ifndef MANDMIMPINT
			if(callType==UPDATE_MOMENTUM_CALL)
			{	// add force (if any) to momentum change
				AddScaledVector(&delPi, &fImp, timestep);
				
				// Add interface energy. (Legacy units g-mm^2/sec^2 or multiply by 1e-9 to get J - kg-m^2/sec^2)
				NodalPoint::interfaceEnergy += rawEnergy;
			}
#else
			if(callType==MASS_MOMENTUM_CALL)
			{	// add total force (in g mm/sec^2) to material field
				mvf[i]->AddFtot(&fImp);
				
				// Add interface energy. (Legacy units g-mm^2/sec^2 or multiply by 1e-9 to get J - kg-m^2/sec^2)
				NodalPoint::interfaceEnergy += rawEnergy;
			}
#endif
		}
		
		// ----------------------------------------------------------------------------------
		// 6. Apply momentum change to the velocity fields
		
		// do not change nodes with boundary conditions (used to only do when postUpdate is true)
		ndptr->AdjustDelPiForBCs(&delPi);
		
		// change momenta
		mvf[i]->ChangeMatMomentum(&delPi,postUpdate,deltime);
		
		// store contact force in rigid particle ftot
		// There are three times this is called in an MPM Step, but to get correct force,
        //    this calculation must only be done during momentum update
        if(callType==UPDATE_MOMENTUM_CALL)
		{	mvf[rigidFld]->AddContactForce(&delPi);
           
            // if conduction on (and here means in momentum update) add frictional heating
            if(mredDE>0. && ConductionTask::matContactHeating)
            {   double qrate = mredDE/massi;
				NodalPoint::frictionWork += qrate;
				conduction->AddFluxCondition(ndptr,fabs(qrate/deltime),true);
            }
        }
	}
}

// Get interfacial contact area and also find tangDel and delt (only used in interface laws)
// Input is ndoptr, voli, volb, delta, norm
// Ouput is tangDel, delt, hperp (optional)
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
	if(hperp!=NULL) *hperp = dist.x;
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
	{	if(MatVelocityField::ActiveField(mvf[i]))
			mvf[i]->CalcVelocityForStrainUpdate();
	}
}

#pragma mark BOUNDARY CONDITIONS

// zero one component of moment and velocity
void CrackVelocityFieldMulti::SetMomVel(Vector *norm)
{	for(int i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
        {	mvf[i]->SetMomentVelocityDirection(norm);
        }
    }
}

// add one component momentum and velocity from BCs
void CrackVelocityFieldMulti::AddMomVel(Vector *norm,double vel)
{	int i;
    for(i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
        {	mvf[i]->AddMomentVelocityDirection(norm,vel);
        }
    }
}

// Reflect one component of velocity and momentum from a node
void CrackVelocityFieldMulti::ReflectMomVel(Vector *norm,CrackVelocityField *rcvf,double reflectRatio)
{	int i;
	MatVelocityField **rmvf = rcvf->GetMaterialVelocityFields();
    for(i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	if(MatVelocityField::ActiveNonrigidField(rmvf[i]))
			{	double rvel = -reflectRatio*DotVectors(norm,&rmvf[i]->pk)/rmvf[i]->mass;
				mvf[i]->AddMomentVelocityDirection(norm,rvel);
			}
		}
	}
}

// set force in direction norm to -p(interpolated)/time such that updated momentum
//    of pk.i + deltime*ftot.i will be zero along norm
void CrackVelocityFieldMulti::SetFtotDirection(Vector *norm,double deltime,Vector *freaction)
{	int i;
    for(i=0;i<maxMaterialFields;i++)
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
void CrackVelocityFieldMulti::ReflectFtotDirection(Vector *norm,double deltime,CrackVelocityField *rcvf,double reflectRatio,Vector *freaction)
{	MatVelocityField **rmvf = rcvf->GetMaterialVelocityFields();
    for(int i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	if(MatVelocityField::ActiveNonrigidField(rmvf[i]))
			{	double rvel = -reflectRatio*DotVectors(norm,&rmvf[i]->pk)/rmvf[i]->mass;
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
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			mass += mvf[i]->mass;
	}
	return mass;
}

// Count number of materials in each crack velocity field on this node
// Sum mass (only of non-rigid materials) and return the result
// Note: for mass, this does not count materials that ignore cracks unless field [0].
// Copy momenta to be restored after force extrapolation
// Uses: Only called once per times step in post extrapolation phase of mass an momentum task
double CrackVelocityFieldMulti::GetTotalMassAndCount(void)
{	
	double mass = 0.;
	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	numberMaterials++;
			if(!mvf[i]->IsRigidField())
			{	mass += mvf[i]->mass;
				
				// copy the extrapolated momenta
				mvf[i]->xpic[PK_COPY] = mvf[i]->pk;
			}
		}
	}
	return mass;
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
	{	if(MatVelocityField:: ActiveNonrigidField(mvf[i]))
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

// get center of mass momentum for all nonrigid material fields in this crack velocity field
// Only counts materials that account for cracks
Vector CrackVelocityFieldMulti::GetCMatMomentum(bool &hasParticles,double *foundMass) const
{	Vector pk;
	ZeroVector(&pk);
	hasParticles = false;
	*foundMass = 0;
	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	AddVector(&pk,&mvf[i]->pk);
			*foundMass += mvf[i]->mass;
			hasParticles = true;
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
	{	if(MatVelocityField:: ActiveNonrigidField(mvf[i]))
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

// get first active rigid field or return NULL. Also return number in rigidFieldNum
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

/* in response to crack contact, change momentum by changing velocity of all 
	nonrigid materials the same amount (and only materials that account
	for cracks)
 
   Change velocity by dP/M, where M is total mass
   Material i velocity becomes vi = pi/mi + dP/M
   Material i momentum change is mi vi = pi + mi dP/M
*/
void CrackVelocityFieldMulti::ChangeCrackMomentum(Vector *delP,bool postUpdate,double deltime)
{
	int i;
	
	// special case for only one material (and it must be nonrigid)
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
		double totMass = GetTotalMass(true);
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


	
