/********************************************************************************
	CrackVelocityFieldMulti.cpp
	nairn-mpm-fea

	Created by John Nairn on 21 August 2009.
	Copyright (c) 2009 John A. Nairn, All rights reserved.
********************************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Nodes/CrackVelocityFieldMulti.hpp"
#include "Nodes/MaterialInterfaceNode.hpp"
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
// throws CommonException() on memory error
void CrackVelocityFieldMulti::MatchMatVelocityFields(MatVelocityField **rmvf)
{	for(int i=0;i<maxMaterialFields;i++)
	{	if(rmvf[i]==NULL) continue;
        
		if(mvf[i]==NULL)
		{	mvf[i]=new MatVelocityField(rmvf[i]->GetFlags());
			if(mvf[i]==NULL) throw CommonException("Memory error allocating material velocity field.",
											   "CrackVelocityFieldMulti::MatchMatVelocityFields");
		}
		
		mvf[i]->Zero();
	}
	numberMaterials=0;
	numberRigidPoints=0;
}

#pragma mark TASK 1 AND 6 METHODS

// Called in intitation to preallocate material velocituy fields
// throws CommonException() on memory error
void CrackVelocityFieldMulti::AddMatVelocityField(int matfld)
{	if(mvf[matfld]==NULL)
	{   mvf[matfld]=new MatVelocityField(MaterialBase::GetMVFFlags(matfld));
        if(mvf[matfld]==NULL)
        {   throw CommonException("Memory error allocating material velocity field.",
												"CrackVelocityFieldMulti::AddMatVelocityField");
        }
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
// Only called by AddTractionForce() and CrackInterfaceForce() when cracks are an imperfect interface
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

	Output changes are only allowed on this node (to be thread safe for parallel)
		changes on mvf[]: pk, vk (if one particle), ftot (if postUpdate is TRUE)
 
	ndptr is parent node to this crack velocity field
	On first call in time step, first and last on pointers to MaterialInterfaceNode * because those
		objects are created for later interface calculations
	postUpdate is TRUE when called between momentum update and particle update and otherwise is false
	throws CommonException() if memory error making interface node
*/
void CrackVelocityFieldMulti::MaterialContactOnCVF(NodalPoint *ndptr,double deltime,int callType,
												   MaterialInterfaceNode **first,MaterialInterfaceNode **last)
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
		RigidMaterialContactOnCVF(rigidMat,multiRigid,ndptr,deltime,callType,first,last);
		return;
	}
	
	// from here on all materials in contact are non-rigid
	AdjustForSymmetry(ndptr,&Pc,false);
	
	// will probably need center of mass if doing displacement check, if imperfect interface,
	// or if transport tasks are doing contact. Because this is amost always, it is always
	// found
	dispc = GetCMDisplacement(ndptr,false);
	ScaleVector(&dispc,1./Mc);
	
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
		bool hasDisplacements = false;
		double dotn,massi=mvf[i]->mass,massRatio=massi/Mc,mred,deln,delt;
        double voli=mvf[i]->GetContactVolume();
		double contactArea = -1.,mredDE = -1.;
		
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
		{	// first look for conditions to ignore contact and interface at this node
			// not sure if these should be ignored for transport contact?
            
			// 1. check nodal volume (this is turned off by setting the materialContactVmin to zero)
			//    (warning: 2D must set grid thickness if it is not 1)
			if(totalVolume/mpmgrid.GetCellVolume()<contact.materialContactVmin) break;
		
			// 2. ignore very small mass nodes - may not be needed
			if(massRatio<1.e-6 || massRatio>0.999999) continue;
		
			// 3. go through contact conditions; break if not in contact or
			//    set inContact to true and break if is in contact. Note that
			//    imperfect interfaces will always proceed to calculations, but it
			//    needs normal vector and needs to know if in contact (to know
			//    whether to use Dnc or Dnt for normal stiffness)
			while(true)
			{	// find -mi(vi-vc) = (ma/mc)pc-pi or momentum change to match ctr of mass momentum
                CopyScaleVector(&delPi,&mvf[i]->pk,-1.);
				AdjustForSymmetry(ndptr, &delPi, false);
                AddScaledVector(&delPi,&Pc,massRatio);

                // Get normal vector by various options
                switch(contact.materialNormalMethod)
                {	case MAXIMUM_VOLUME_GRADIENT:
                    {	// Use mat with largest magnitude volume gradient
                        Vector normi;
                        GetVolumeGradient(i,ndptr,&normi,1.);
                        
                        // compare magnitude of volume gradients
                        double magi=sqrt(DotVectors(&normi,&normi));
                        double magj=sqrt(DotVectors(&otherGrad,&otherGrad));
                        if(magi >= magj)
                            CopyScaleVector(&norm,&normi,1./magi);		// use material i
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

                // get approach direction momentum from delPi.n (actual (vc-vi) = delPi/mi)
                dotn=DotVectors(&delPi,&norm);
                
                // 3a. With this check, any movement apart will be taken as noncontact
                // Also, frictional contact assumes dotn<0
                if(dotn>=0.) break;
                
                // 3b. Displacement check: Get COD info or delta = dispcScaled - dispi
                if(contact.displacementCheck)
				{	// 4. get other mass and ignore if very small mass in other materials
                    mred=(Mc-massi)/Mc;
                    if(mred<1.e-6) break;
                    
					// get COD
					mred = GetContactCOD(ndptr,&delta,&dispc,&mvf[i]->disp,massi,mred,&hasDisplacements);
				
					// to get normal velocity delta v = (Mc/(Mc-mi)) (delta p/mi)
					double dvel = dotn/mred;
                    
					// 5. check for contact
					if(contact.MaterialContact(&delta,&norm,dvel,postUpdate,deltime)==SEPARATED) break;
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
			
			// if not in contact, and not imperfect interface (which needs to contine), then
			// done if only two materials (and the two use same normal) or continue to next materiasl
			if(theContactLaw->ContactIsDone(inContact))
			{	if(numberMaterials==2 && contact.materialNormalMethod!=EACH_MATERIALS_MASS_GRADIENT) break;
				continue;
			}
		}
		
		// Here inContact is true, but if law is IMPERFECT_INTERFACE, may be true or false
        // adjust momentum change as needed
		if(comContact)
		{	// use delPi found above and no further change needed
		}
		else if(theContactLaw->IsFrictionalContact())
		{	bool getHeating = (callType==UPDATE_MOMENTUM_CALL) && ConductionTask::matContactHeating;
			mred = massi*(Mc-massi)/Mc;
			if(contactArea<0. && theContactLaw->FrictionLawNeedsContactArea())
			{	contactArea = GetContactArea(ndptr,voli,GetVolumeNonrigid(false)-voli,
											 &delta,&norm,&tangDel,&deln,&delt);
			}
			
			// Find (ma Fb - mb Fa)/Mc = (ma Fc - Mc Fa)/Mc for 2nd order heating
			Vector at = GetCMatFtot();
			ScaleVector(&at, massi);
			AddScaledVector(&at, mvf[i]->GetFtotPtr(), -Mc);
			ScaleVector(&at,1./Mc);
			
			if(!theContactLaw->GetFrictionalDeltaMomentum(&delPi,&norm,dotn,&mredDE,mred,getHeating,contactArea,inContact,deltime,&at))
			{	if(numberMaterials==2 && contact.materialNormalMethod!=EACH_MATERIALS_MASS_GRADIENT) break;
				continue;
			}
		}
		else
		{	// may not have found displacement above yet
			if(!hasDisplacements)
			{	mred = (Mc-massi)/Mc;
				mred = GetContactCOD(ndptr,&delta,&dispc,&mvf[i]->disp,massi,mred,&hasDisplacements);
			}
			
			// get interface forces for future use and get non-uniform grid correction
			Vector fImp;
			double rawEnergy;
			if(contactArea<0.)
			{	contactArea = GetContactArea(ndptr,voli,GetVolumeNonrigid(false)-voli,
											 &delta,&norm,&tangDel,&deln,&delt);
			}
			bool createNode = theContactLaw->GetInterfaceForcesForNode(&norm,&fImp,&rawEnergy,
														contactArea,&delPi,dotn,inContact,postUpdate,mred,
														&tangDel,deln,delt);
			
			if(createNode && first!=NULL)
			{	// decide if force balance can get other node too
				int iother = numberMaterials==2 && contact.materialNormalMethod!=EACH_MATERIALS_MASS_GRADIENT ? ipaired : -1 ;
				
				// only add interface energy once (i.e. when more than 2 or using own grad)
				if(iother==-1 && hasInterfaceEnergy) rawEnergy = 0.;
				
				// create node to add internal force later, if needed set first one
				*last = new MaterialInterfaceNode(ndptr,fieldNum,i,iother,&fImp,rawEnergy,*last);
				if(*last == NULL)
				{	throw CommonException("Memory error allocating storage for a material interface node.",
										  "CrackVelocityFieldMulti::MaterialContactOnCVF");
				}
				if(*first==NULL) *first=*last;
				
				// has energy at least once
				hasInterfaceEnergy = true;
			}
		}
		
		// on post update contact, do not change nodes with boundary conditions
		unsigned char fixedDirection=ndptr->fixedDirection;
        if(postUpdate && (fixedDirection&XYZ_SKEWED_DIRECTION))
		{	if(fixedDirection&X_DIRECTION) delPi.x=0.;
			if(fixedDirection&Y_DIRECTION) delPi.y=0.;
			if(fixedDirection&Z_DIRECTION) delPi.z=0.;
		}
		
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
// means exactly one is a rigid material
//	(no rigid materials handled in MaterialContactOnCVF(), two rigid materials is an error)
// throws CommonException() if memory error making interface node
void CrackVelocityFieldMulti::RigidMaterialContactOnCVF(int rigidFld,bool multiRigid,NodalPoint *ndptr,double deltime,int callType,
												   MaterialInterfaceNode **first,MaterialInterfaceNode **last)
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
	if(mvf[rigidFld]->numberPoints==1)
        rvel = mvf[rigidFld]->GetVelocity();
	else
		CopyScaleVector(&rvel, &mvf[rigidFld]->pk, 1./actualRigidVolume);
	AdjustForSymmetry(ndptr,&rvel,false);
	
	// loop over each material (skipping the one rigid material)
	bool postUpdate = callType != MASS_MOMENTUM_CALL;
	for(i=0;i<maxMaterialFields;i++)
    {	if(!MatVelocityField::ActiveField(mvf[i]) || mvf[i]->IsRigidField()) continue;
        
		// First determine contact law with rigid material
		ContactLaw *theContactLaw = contact.GetMaterialContactLaw(i,rigidFld);
        
        // get contact parameters
        Vector delta,tangDel;
        bool inContact = false;
        bool hasDisplacements = false;
	
		// some variables
		Vector norm,delPi;
		double dotn,massi=mvf[i]->mass,voli=mvf[i]->GetContactVolume(),deln,delt;
		double contactArea = -1.,mredDE = -1.;
        
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
        {   // first look for conditions to ignore contact and interface at this node
            
            // 1. check nodal volume (this is turned off by setting the materialContactVmin to zero)
            //    (warning: 2D must set grid thickness if it is not 1)
            if(GetVolumeTotal(ndptr)/mpmgrid.GetCellVolume()<contact.materialContactVmin) continue;
            
            // 2. ignore very small interactions
            double volRatio=voli/(rigidVolume+GetVolumeNonrigid(false));
            if(volRatio<.001 || volRatio>0.999) continue;
            //if(volRatio<1.e-6 || volRatio>0.999999) continue;
            
            // 3. go through contact conditions; break if not in contact or
            //    set inContact to true and break if is in contact. Note that
            //    imperfect interfaces will always proceed to calculations, but it
            //    needs normal vector and needs to know if in contact (to know
            //    whether to use Dnc or Dnt for normal stiffness)
            while(true)
			{	// Get normal vector by various options
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

                 // get approach direction momentum form delPi.n (actual (vr-vi).n = delPi.n/mi)
                dotn=DotVectors(&delPi,&norm);
            
                // 3a. With this check, any movement apart will be taken as noncontact
                // Also, frictional contact assumes dotn<0
                if(dotn>=0.) break;
            
                // 3b. Displacement check
                if(contact.displacementCheck)
                {	// rigid material displacement was scaled by volume, while non-rigid was weighted by mass
					Vector dispi;
					CopyScaleVector(&delta,&mvf[rigidFld]->disp,1./actualRigidVolume);
					AdjustForSymmetry(ndptr,&delta,false);
                    CopyScaleVector(&dispi,&mvf[i]->disp,1./massi);
					AdjustForSymmetry(ndptr,&dispi,false);
					SubVector(&delta,&dispi);
					hasDisplacements = true;
                
                    // convert dotn to velocity of approach
                    double dvel = dotn/massi;
                
                    // check for contact
                    if(contact.MaterialContact(&delta,&norm,dvel,postUpdate,deltime)==SEPARATED) break;
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
			
            // continue if not in contact, unless it is an imperfect interface law
            if(theContactLaw->ContactIsDone(inContact)) continue;
        }
		
		if(comContact)
		{	// use delPi found above and no further change needed
		}
		else if(theContactLaw->IsFrictionalContact())
		{	bool getHeating = (callType==UPDATE_MOMENTUM_CALL) && ConductionTask::matContactHeating;
			if(contactArea<0. && theContactLaw->FrictionLawNeedsContactArea())
			{	contactArea = GetContactArea(ndptr,voli,rigidVolume,&delta,&norm,&tangDel,&deln,&delt);
			}
			
			// Find -fi,a for second order heating
			Vector at = mvf[i]->GetFtot();
			ScaleVector(&at, -1.);
			
			if(!theContactLaw->GetFrictionalDeltaMomentum(&delPi,&norm,dotn,&mredDE,massi,getHeating,contactArea,inContact,deltime,&at))
				continue;
		}
		else
		{	// get displacement of material i and the rigid material
			if(!hasDisplacements)
			{   // rigid material displacement was scaled by volume, while non-rigid was weighted by mass
				Vector dispi;
				CopyScaleVector(&delta,&mvf[rigidFld]->disp,1./actualRigidVolume);
				AdjustForSymmetry(ndptr,&delta,false);
				CopyScaleVector(&dispi,&mvf[i]->disp,1./massi);
				AdjustForSymmetry(ndptr,&dispi,false);
				SubVector(&delta,&dispi);
				hasDisplacements = true;
			}
			
			// get interface forces for future use and nonuniform grid correction
			Vector fImp;
			double rawEnergy;
			if(contactArea<0.)
			{	contactArea = GetContactArea(ndptr,voli,rigidVolume,&delta,&norm,&tangDel,&deln,&delt);
			}
			bool createNode = theContactLaw->GetInterfaceForcesForNode(&norm,&fImp,&rawEnergy,
														contactArea,&delPi,dotn,inContact,postUpdate,massi,
														&tangDel,deln,delt);
			
			if(createNode && first!=NULL)
			{	// create node to add internal force later, if needed set first one
				*last = new MaterialInterfaceNode(ndptr,fieldNum,i,-1,&fImp,rawEnergy,*last);
				if(*last==NULL)
				{	throw CommonException("Memory error allocating storage for a material interface node.",
										  "CrackVelocityFieldMulti::RigidMaterialContactOnCVF");
				}
				if(*first==NULL) *first = *last;
			}
		}
		
		// on post update contact, do not change nodes with boundary conditions
		unsigned char fixedDirection=ndptr->fixedDirection;
		if(postUpdate && (fixedDirection&XYZ_SKEWED_DIRECTION))
		{	if(fixedDirection&X_DIRECTION) delPi.x=0.;
			if(fixedDirection&Y_DIRECTION) delPi.y=0.;
			if(fixedDirection&Z_DIRECTION) delPi.z=0.;
		}
		
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

// Get interfacial contact area and also find tangDel, deln, and delt
double CrackVelocityFieldMulti::GetContactArea(NodalPoint *ndptr,double voli,double volb,Vector *delta,
											   Vector *norm,Vector *tangDel,double *delnout,double *deltout) const
{
	// Get raw surface area, it is divided by hperp in GetInterfaceForceForNode()
	// Scale voltot=voli+volb to voltot*sqrt(2*vmin/voltot) = sqrt(2*vmin*vtot)
	// multiple by position for axisymmetric
	double rawSurfaceArea = sqrt(2.*fmin(voli,volb)*(voli+volb));
	if(fmobj->IsAxisymmetric()) rawSurfaceArea *= ndptr->x;			// times position if axisym
	
	// perpendicular distance to correct contact area and contact by positions
    double dist;
    
    // normal displacement (norm is normalized) = delta . norm, subtract adjustment when using position
	// which have been precalculated
    double deln = DotVectors(delta,norm);
	
    // tangential vector in tangDel
    CopyVector(tangDel,delta);
    AddScaledVector(tangDel,norm,-deln);				// delta - deln (n) = dott (t)
    double delt=sqrt(DotVectors(tangDel,tangDel));
    if(!DbleEqual(delt,0.)) ScaleVector(tangDel,1/delt);
	
	// if using displacements, find dist  with tangDel
	if(mpmgrid.GetContactByDisplacements())
	{	dist = mpmgrid.GetPerpendicularDistance(norm, tangDel, delt);
	}
	
	else
	{	// Cannot use tangDEl here, which means 3D non regular will be less accurate
		// But, all regular and all 2D will be correct without tang
		// Future Goal: get hperp without needing tang for general 3D case
		// (for efficiency, call hperp method separately)
		dist = mpmgrid.GetPerpendicularDistance(norm, NULL, 0.);
        deln -= mpmgrid.positionCutoff*dist;
	}
 	
	*delnout = deln;
	*deltout = delt;
	return rawSurfaceArea/dist;
}

// scale displacements to get delta = (Mc/(Mc-mi))*dispc - (Mc/(Mc-mi))*(mvf[i]->disp/mi)
// dispi is displacement (or position) for material i and (Mc/(Mc-mi))*dispc is displacement
//    (or position) for the virtual paired material (combination of other materials)
// the separation vector is their difference
// return reduced mass = (mi*(Mc-mi))/Mc while on call mred = (Mc-mi)/Mc
double CrackVelocityFieldMulti::GetContactCOD(NodalPoint *ndptr,Vector *delta,Vector *dispc,Vector *mvfdisp,
											double mi,double mred,bool *hasDisplacements)
{
	Vector dispi;
	CopyScaleVector(delta,dispc,1./mred);
	mred *= mi;                // now mred = (mi*(Mc-mi))/Mc
	CopyScaleVector(&dispi,mvfdisp,1./mred);
	AdjustForSymmetry(ndptr,&dispi,false);
	SubVector(delta,&dispi);
	*hasDisplacements=true;
	return mred;
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
void CrackVelocityFieldMulti::ReflectMomVel(Vector *norm,CrackVelocityField *rcvf)
{	int i;
	MatVelocityField **rmvf = rcvf->GetMaterialVelocityFields();
    for(i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	if(MatVelocityField::ActiveNonrigidField(rmvf[i]))
			{	double rvel = -DotVectors(norm,&rmvf[i]->pk)/rmvf[i]->mass;
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
void CrackVelocityFieldMulti::ReflectFtotDirection(Vector *norm,double deltime,CrackVelocityField *rcvf,Vector *freaction)
{	MatVelocityField **rmvf = rcvf->GetMaterialVelocityFields();
    for(int i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
		{	if(MatVelocityField::ActiveNonrigidField(rmvf[i]))
			{	double rvel = -DotVectors(norm,&rmvf[i]->pk)/rmvf[i]->mass;
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
// Uses: Only called once per times step in post extrapolation phase of mass an momentum task
double CrackVelocityFieldMulti::GetTotalMassAndCount(void)
{	
	double mass = 0.;
	for(int i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	numberMaterials++;
			if(!mvf[i]->IsRigidField())
				mass += mvf[i]->mass;
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

// add contact force on rigid material to the input vector
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

/* in response to crack contact, change moment by changing velocity of all 
	nonrigid materials the same amount
 
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


	
