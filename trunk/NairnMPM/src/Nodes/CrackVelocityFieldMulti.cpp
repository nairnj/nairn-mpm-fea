/********************************************************************************
	CrackVelocityFieldMulti.cpp
	NairnMPM

	Created by John Nairn on 21 August 2009.
	Copyright (c) 2009 John A. Nairn, All rights reserved.
********************************************************************************/

#include "NairnMPM_Class/NairnMPM.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Nodes/CrackVelocityFieldMulti.hpp"
#include "Nodes/MaterialInterfaceNode.hpp"
#include "Exceptions/MPMTermination.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"

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

// add "mass" for rigid particle (task 1) - it counts particles too
void CrackVelocityFieldMulti::AddMassTask1(int matfld,double mnode)
{	mvf[matfld]->mass += mnode;
	numberRigidPoints++;
}

// Add to mass gradient
// This gradient is only used in multimaterial contact calculations
// It may be wrong for axisymmetric, but that is fixed in GetVolumeGradient(). It is
//		faster there because that is called less while this is called for every node-particle pair
void CrackVelocityFieldMulti::AddVolumeGradient(int matfld,MPMBase *mptr,double dNdx,double dNdy,double dNdz)
{	double Vp = mptr->GetVolume(DEFORMED_AREA);
	mvf[matfld]->volumeGrad->x+=Vp*dNdx;
	mvf[matfld]->volumeGrad->y+=Vp*dNdy;
	mvf[matfld]->volumeGrad->z+=Vp*dNdz;
}

// Count number of materials, and return total mass
// Get mass only of non-rigid materials, but count them all
double CrackVelocityFieldMulti::GetTotalMassAndCount(void)
{	int i;
	double mass=0.;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
		{	numberMaterials++;
			if(!mvf[i]->rigidField)
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
	
	// add unscaled volume to this crack velocity field (deprecated, need to fix)
	//unscaledRigidVolume+=cvfm->UnscaledVolumeRigid();
	
	// sum momentum, displacement, and mass grad (velocity is same) into material velocity field
	mvf[rigidFieldNum]->numberPoints+=rmvf->numberPoints;
	AddVector(&mvf[rigidFieldNum]->pk,&rmvf->pk);
	AddVector(&mvf[rigidFieldNum]->disp,&rmvf->disp);
	AddVector(mvf[rigidFieldNum]->volumeGrad,rmvf->volumeGrad);
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
	unscaledRigidVolume-=UnscaledVolumeRigid();
	*/
	
	// create material field if needed
	int initialRigidPoints=0;
	//double initialRigidVolume=0.;
	if(mvf[rigidFieldNum]==NULL)
	{	mvf[rigidFieldNum]=new MatVelocityField(TRUE);
		if(mvf[rigidFieldNum]==NULL) throw CommonException("Memory error allocating material velocity field.",
													"CrackVelocityFieldMulti::CopyRigidFrom");
		numberMaterials++;					// just added a material to this crack velocity field
	}
	else
	{	initialRigidPoints=mvf[rigidFieldNum]->numberPoints;
		// initialRigidVolume=UnscaledVolumeRigid(); deprecated need to fix
		if(initialRigidPoints==0) numberMaterials++;
	}
	
	// reference to source field
	MatVelocityField *rmvf = cvfm->mvf[rigidFieldNum];
	
	// add number rigid points this crack velocity field
	numberRigidPoints += rmvf->numberPoints-initialRigidPoints;
	numberPoints += rmvf->numberPoints-initialRigidPoints;
	
	// add unscaled volume to this crack velocity field (may be wrong due to recent change in unscaled volumes)
	//unscaledRigidVolume+=cvfm->UnscaledVolumeRigid()-initialRigidVolume; deprecated need to fix
	
	// copy momentum, displacement, and mass grad (velocity is same) into material velocity field
	mvf[rigidFieldNum]->numberPoints = rmvf->numberPoints;
	CopyVector(&mvf[rigidFieldNum]->pk,&rmvf->pk);
    mvf[rigidFieldNum]->SetVelocity(rmvf->GetVelocityPtr());
	CopyVector(&mvf[rigidFieldNum]->disp,&rmvf->disp);
	CopyVector(mvf[rigidFieldNum]->volumeGrad,rmvf->volumeGrad);
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
			{	mvf[i]->AddFint(f);
				break;
			}
		}
	}
	
	// more than one material, add to nonrigid materials only
	else
	{	double totMass=GetTotalMass();
		for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
				mvf[i]->AddFint(f,mvf[i]->mass/totMass);
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
			{	mvf[i]->AddFext(f);
				break;
			}
		}
	}
	
	// more than one material
	else
	{	double totMass=GetTotalMass();
		for(i=0;i<maxMaterialFields;i++)
		{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
				mvf[i]->AddFext(f,mvf[i]->mass/totMass);
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
            mvf[i]->UpdateMomentum(timestep);
    }
}

#pragma mark TASK 6 METHODS

// zero momentum and displacement at a node for new calculations
// but can do the calculation for rigid particles here
void CrackVelocityFieldMulti::RezeroNodeTask6(double deltaTime)
{	int i;
    for(i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveField(mvf[i]))
		{	if(!mvf[i]->rigidField)
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

#pragma mark MATERIAL CONTACT

/* Called in multimaterial mode to check contact at nodes with multiple materials

	Input parameters:
		mvf[]->mass,pk,volumeGrad,disp (if contact by displacements)

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
		{	if(!mvf[i]->rigidField)
			{	// real material
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
	bool hasInterfaceEnergy = FALSE;
	for(i=0;i<maxMaterialFields;i++)
    {	if(!MatVelocityField::ActiveField(mvf[i])) continue;
		
		// some variables
		Vector norm,delPi,dispcScaled,dispi;
        bool hasDisplacements = FALSE;
		double dotn,massi=mvf[i]->mass,massRatio=massi/Mc,mred;
        double voli=mvf[i]->GetContactVolume();
		
		// First determine contact law from other material with most volume
        // and find total other volume and volume weighted mean volume gradient
		double maxOtherMaterialVolume=0.;
		int ipaired=0;
		for(j=0;j<maxMaterialFields;j++)
		{	if(j==i || !MatVelocityField::ActiveField(mvf[j])) continue;
			double matVolume=mvf[j]->GetContactVolume();
			if(matVolume>maxOtherMaterialVolume)
			{	maxOtherMaterialVolume=matVolume;
				ipaired=j;
			}
		}
		
		// problem if ipaired not found, but it will be found
		int maxContactLaw=contact.GetMaterialContactLaw(i,ipaired);
		double maxFriction=0.,Dn,Dnc,Dt;
        bool inContact = FALSE;
        
        // get contact parameters
        if(maxContactLaw == IMPERFECT_INTERFACE)
        {   // fetch interface properties
            contact.GetMaterialInterface(i,ipaired,&Dn,&Dnc,&Dt);
        }
        else
            maxFriction=contact.GetMaterialFriction(i,ipaired);
		
        // If needed, determine if the surfaces are in contact
        if(maxContactLaw == NOCONTACT)
        {   // find -mi(vi-vc) = (ma/mc)pc-pi or momentum change to match ctr of mass momentum
            CopyScaleVector(&delPi,&mvf[i]->pk,-1.);
            AddScaledVector(&delPi,&Pc,massRatio);
            inContact = TRUE;
        }
        
        else
        {   // first look for conditions to ignore contact and interface at this node
            
            // 1. check nodal volume (this is turned off by setting the materialContactVmin to zero)
            //    (warning: 2D must set grid thickness if it is not 1)
            if(GetVolumeTotal(nd[nodenum]->x)/mpmgrid.GetCellVolume()<contact.materialContactVmin) continue;
            
            // 2. ignore very small mass nodes - may not be needed
            if(massRatio<1.e-6 || massRatio>0.999999) continue;
            
            // 3. go through contact conditions; break if not in contact or
            //    set inContact to true and break if is in contact. Note that
            //    imperfect interfaces will always proceed to calculations, but it
            //    needs normal vector and needs to know if in contact (to know
            //    whether to use Dnc or Dnt for normal stiffness)
            while(TRUE)
            {   // find -mi(vi-vc) = (ma/mc)pc-pi or momentum change to match ctr of mass momentum
                CopyScaleVector(&delPi,&mvf[i]->pk,-1.);
                AddScaledVector(&delPi,&Pc,massRatio);

                // Get normal vector by various options
                switch(contact.materialNormalMethod)
                {	case MAXIMUM_VOLUME_GRADIENT:
                    {	// Use mat with largest magnitude volume gradient
                        Vector normi,normj;
                        nd[nodenum]->GetVolumeGradient(vfld,i,&normi,1.);
                        nd[nodenum]->GetVolumeGradient(vfld,ipaired,&normj,-1.);
                        
                        // compare magnitude of volume gradients
                        double magi=sqrt(DotVectors(&normi,&normi));
                        double magj=sqrt(DotVectors(&normj,&normj));
                        if(magi >= magj)
                            CopyScaleVector(&norm,&normi,1./magi);		// use material i
                        else
                            CopyScaleVector(&norm,&normj,1./magj);		// use material j
                        break;
                    }
                        
                    case MAXIMUM_VOLUME:
						// Use mat with most volume
                        if(voli >= maxOtherMaterialVolume)
                            nd[nodenum]->GetVolumeGradient(vfld,i,&norm,1.);
                        else
                            nd[nodenum]->GetVolumeGradient(vfld,ipaired,&norm,-1.);
                        CopyScaleVector(&norm,&norm,1./sqrt(DotVectors(&norm,&norm)));
                       break;
                        
                    case AVERAGE_MAT_VOLUME_GRADIENTS:
                    {	// get mass gradients material i and for other material(s)
                        Vector normi,normj,otherGrad;
                        nd[nodenum]->GetVolumeGradient(vfld,i,&normi,1.);
                        
                        // Find total other volume weighted volume gradient
                        //nd[nodenum]->GetVolumeGradient(vfld,ipaired,&normj,-1.);        // to just use paired one
                        ZeroVector(&otherGrad);
                        for(j=0;j<maxMaterialFields;j++)
                        {	if(j==i || !MatVelocityField::ActiveField(mvf[j])) continue;
                             
                            // Finding -V grad V = -Sum V_j grad V_j 
                            nd[nodenum]->GetVolumeGradient(vfld,j,&normj,-1.);
                            AddScaledVector(&otherGrad,&normj,mvf[j]->GetContactVolume());
                        }
                        // No need to divide by otherVolume to normalize because we want
                        // avg grad * other volume, which was found above
                       
                        // volume weighted mean of volume gradients
                        //  = (voli * grad voli + otherVolume * grad otherVolume)/(total volume)
                        CopyScaleVector(&norm,&normi,voli);
                        AddScaledVector(&norm,&otherGrad,1.);
                        
                        // normalize
                        double magi=sqrt(DotVectors(&norm,&norm));
                        ScaleVector(&norm,1./magi);
                        break;
                    }
                        
                    case EACH_MATERIALS_MASS_GRADIENT:
                        // Use each material's own gradient and handle separately (i.e. does not conserve momentum)
                        nd[nodenum]->GetVolumeGradient(vfld,i,&norm,1.);
                        ScaleVector(&norm,1./sqrt(DotVectors(&norm,&norm)));
                        break;
                        
                   default:
                        break;
                }
                
                // current options: 3 (give an axis of z rotation angle), 5 (spherical particle, only here in nonrigid)
 				if(fmobj->dflag[0]==3)
				{   // normal along +/-x, +/-y or +/-z from flag[1] = +/-1, +/-2, or +/-3
                    //    otherwisse rotates cw about z axis by that number of degrees n = (cos(angle),-sin(angle),0)
					// This should be the normal vector pointing out of lower numbered material
					int normAxis = fmobj->dflag[1];
					if(normAxis==1 || normAxis==-1)
					{   norm.x = (double)normAxis;
						norm.y = 0.;
						norm.z = 0.;
					}
					else if(normAxis==2 || normAxis==-2)
					{   norm.x = 0.;
						norm.y = normAxis>0 ? 1. : -1. ;
						norm.z = 0.;
					}
					else if(normAxis==3 || normAxis==-3)
					{   norm.x = 0.;
						norm.y = 0.;
						norm.z = normAxis>0 ? 1. : -1. ;
					}
                    else
                    {   double radAngle=(double)fmobj->dflag[1]*PI_CONSTANT/180.;
                        norm.x=cos(radAngle);
                        norm.y=-sin(radAngle);
                        norm.z = 0.;
                    }
				}
				else if(fmobj->dflag[0]==5)
				{	norm.x = nd[nodenum]->x;
					norm.y = nd[nodenum]->y;
					norm.z = nd[nodenum]->z;
					double normmag=sqrt(DotVectors(&norm,&norm));
					norm.x /= normmag;
					norm.y /= normmag;
					norm.z /= normmag;
				}
                
                // get approach direction momentum from delPi.n (actual (vc-vi) = delPi/mi)
                dotn=DotVectors(&delPi,&norm);
                
                // 3a. With this check, any movement apart will be taken as noncontact
                // Also, frictional contact assumes dotn<0
                if(dotn>=0.) break;
                
                // 3b. Displacement check
                if(contact.displacementCheck)
                {	// 4. get other mass and ignore if very small mass in other materials
                    mred=(Mc-massi)/Mc;
                    if(mred<1.e-6) break;
                    
                    // scale displacements to get delta = (Mc/(Mc-mi))*disp - (Mc/(Mc-mi))*(mvf[i]->disp/mi)
                    // dispi is displacement (or position) for material i and dispcScaled is displacement
                    //    (or position) for the virtual paired material (combination of other materials)
                    // the separation vector is their difference
                    CopyScaleVector(&dispcScaled,&dispc,1./mred);
                    mred*=massi;                // now (ma mother)/(ma+mother)
                    CopyScaleVector(&dispi,&mvf[i]->disp,1./mred);
                    hasDisplacements = TRUE;
                    
                    // to get normal velocity delta v = (Mc/(Mc-mi)) (delta p/mi)
                    double dvel = dotn/mred;
                    
                    // 5. check for contact
                    if(contact.MaterialContact(&dispi,&dispcScaled,&norm,dvel,postUpdate,deltime)==SEPARATED) break;
                }
                
                // passed all tests
                inContact = TRUE;
                break;
            }
            
            // continue if not in contact, unless it is an imperfect interface law
            if(!inContact && (maxContactLaw!=IMPERFECT_INTERFACE))
			{	// special case two materials for efficiency (and if both will find normal the same way)
				// can break out of contact because now know to be no contact, otherwise continue to next material
				if(numberMaterials==2 && contact.materialNormalMethod!=EACH_MATERIALS_MASS_GRADIENT)
					break;
				continue;
			}
		}
		
        // adjust momentum change as needed
		switch(maxContactLaw)
		{	case STICK:
			case NOCONTACT:
                // use delPi as is to get center of mass motion
                // NOCONTACT is here always
                // STICK here only if in contact
				break;
				
			case FRICTIONLESS:
				CopyScaleVector(&delPi,&norm,dotn);
				break;
				
			case FRICTIONAL:
                GetFrictionalDeltaMomentum(&delPi,&norm,dotn,maxFriction);
				break;
            
            case IMPERFECT_INTERFACE:
			{	// get displacement of material i and the virtual opposite material
				if(!hasDisplacements)
				{   if(!contact.displacementCheck)
					{	dispc=GetCMDisplacement();
						ScaleVector(&dispc,1./Mc);
					}
					
					// scale displacements to get delta = (Mc/(Mc-mi))*dispc - (Mc/(Mc-mi))*(mvf[i]->disp/mi)
					//    of delta = disp(virtual) - dist(mat i) = dispcScaled - dispi
					mred = (Mc-massi)/Mc;
					CopyScaleVector(&dispcScaled,&dispc,1./mred);
					mred *= massi;                  // now (ma mother)/(ma+mother)
					CopyScaleVector(&dispi,&mvf[i]->disp,1./mred);
				}
				
				// find displacement difference vector
				Vector delta=dispcScaled;
				SubVector(&delta,&dispi);
				
				// get interface forces for future use and get non-uniform grid correction
				Vector fImp;
				double rawEnergy;
				
				// Get raw surface area, it is divided by hperp in GetInterfaceForceForNode()
				// Scale voltot=voli+volb to voltot*sqrt(2*vmin/voltot) = sqrt(2*vmin*vtot)
				double volb = GetVolumeNonrigid()-voli;
				double rawSurfaceArea = sqrt(2.*fmin(voli,volb)*(voli+volb));
				//double rawSurfaceArea = 2.*fmin(voli,volb);
				//double rawSurfaceArea = (voli+volb);
				
				// multiple by position for axisymmetric
				if(fmobj->IsAxisymmetric()) rawSurfaceArea *= nd[nodenum]->x;
				
				// get forces
				bool createNode = GetInterfaceForcesForNode(&delta,&norm,Dn,Dnc,Dt,&fImp,&rawEnergy,
															rawSurfaceArea,&delPi,dotn,inContact,postUpdate,mred);
				
				if(createNode && !postUpdate)
				{	// decide if force balance can get other node too
					int iother = numberMaterials==2 && contact.materialNormalMethod!=EACH_MATERIALS_MASS_GRADIENT ? ipaired : -1 ;
					
					// only add interface energy once (i.e. when more than 2 or using own grad)
					if(iother==-1 && hasInterfaceEnergy) rawEnergy = 0.;
					
					// create node to add internal force later
					MaterialInterfaceNode::currentNode=new MaterialInterfaceNode(nd[nodenum],vfld,i,iother,&fImp,rawEnergy);
					if(MaterialInterfaceNode::currentNode==NULL)
					{	throw CommonException("Memory error allocating storage for a material interface node.",
													"CrackVelocityFieldMulti::MaterialContact");
					}
					
					// has energy at least once
					hasInterfaceEnergy = TRUE;
				}
                break;
			}
				
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
		if(numberMaterials==2 && contact.materialNormalMethod!=EACH_MATERIALS_MASS_GRADIENT)
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
	// get rigid material volume for contact and actual volume
	// these are the same except in axisymmtric where first is area and second is volume per radian
	double rigidVolume = mvf[rigidFld]->GetContactVolume();
	double actualRigidVolume = mvf[rigidFld]->mass;
	
	// get rigid volume
	// rigid material with position-dependent velocities may have mixed velocities
	Vector rvel;
	if(mvf[rigidFld]->numberPoints==1)
        rvel = mvf[rigidFld]->GetVelocity();
	else
		CopyScaleVector(&rvel, &mvf[rigidFld]->pk, 1./actualRigidVolume);
	
	// loop over each material (skipping the one rigid material)
	int i;
	for(i=0;i<maxMaterialFields;i++)
    {	if(!MatVelocityField::ActiveField(mvf[i]) || i==rigidFld) continue;
		
		// First determine contact law with rigid material
		int maxContactLaw=contact.GetMaterialContactLaw(i,rigidFld);
        
        // get contact parameters
        Vector rigidDisp,dispi;
        double maxFriction,Dn,Dnc,Dt;
        bool inContact = FALSE;
        bool hasDisplacements = FALSE;
        if(maxContactLaw == IMPERFECT_INTERFACE)
        {   // fetch interface properties
            contact.GetMaterialInterface(i,rigidFld,&Dn,&Dnc,&Dt);
        }
        else
            maxFriction=contact.GetMaterialFriction(i,rigidFld);
		
		// some variables
		Vector norm,delPi;
		double dotn,massi=mvf[i]->mass,voli=mvf[i]->GetContactVolume();
        
        // find -mi(vi-vr) = - pi + mi*vr, which is change in momentum to match the rigid particle velocity
        CopyScaleVector(&delPi,&mvf[i]->pk,-1.);
        AddScaledVector(&delPi,&rvel,massi);
        
        // if needed, determine if the surfaces are in contact
        if(maxContactLaw==NOCONTACT)
        {   // will use -mi(vi-vr) = - pi + mi*vr to match rigid particle velocity
            inContact=TRUE;
        }
        
        else
        {   // first look for conditions to ignore contact and interface at this node
            
            // 1. check nodal volume (this is turned off by setting the materialContactVmin to zero)
            //    (warning: 2D must set grid thickness if it is not 1)
            if(GetVolumeTotal(nd[nodenum]->x)/mpmgrid.GetCellVolume()<contact.materialContactVmin) continue;
            
            // 2. ignore very small interactions
            double volRatio=voli/(rigidVolume+GetVolumeNonrigid());
            if(volRatio<.001 || volRatio>0.999) continue;
            //if(volRatio<1.e-6 || volRatio>0.999999) continue;
            
            // 3. go through contact conditions; break if not in contact or
            //    set inContact to true and break if is in contact. Note that
            //    imperfect interfaces will always proceed to calculations, but it
            //    needs normal vector and needs to know if in contact (to know
            //    whether to use Dnc or Dnt for normal stiffness)
            while(TRUE)
			{	// Get normal vector by various options
                switch(contact.materialNormalMethod)
                {	case MAXIMUM_VOLUME_GRADIENT:
                    {	// Use mat with largest magnitude volume gradient
                        Vector normi,normj;
                        nd[nodenum]->GetVolumeGradient(vfld,i,&normi,1.);
                        nd[nodenum]->GetVolumeGradient(vfld,rigidFld,&normj,-1.);
                    
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
                            nd[nodenum]->GetVolumeGradient(vfld,i,&norm,1.);
                        else
                            nd[nodenum]->GetVolumeGradient(vfld,rigidFld,&norm,-1.);
                        CopyScaleVector(&norm,&norm,1./sqrt(DotVectors(&norm,&norm)));
						break;
                        
                    case AVERAGE_MAT_VOLUME_GRADIENTS:
                    {	// get volume-weighted mean of volume gradiants of material and rigid material
                        Vector normi,normj;
                        nd[nodenum]->GetVolumeGradient(vfld,i,&normi,1.);
                        nd[nodenum]->GetVolumeGradient(vfld,rigidFld,&normj,-1.);
                    
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
                        nd[nodenum]->GetVolumeGradient(vfld,i,&norm,1.);
                        CopyScaleVector(&norm,&norm,1./sqrt(DotVectors(&norm,&norm)));
                        break;

                    default:
                        break;
                }			
            
                // current options: 3 (give an axis of z rotation angle)
                //                  4 (cutting, but only here in rigid contact)
                if(fmobj->dflag[0]==3)
                {   // normal along +/-x, +/-y or +/-z from flag[1] = +/-1, +/-2, or +/-3
                    //    otherwisse rotates cw about z axis by that number of degrees n = (cos(angle),-sin(angle),0)
                    // This should be the normal vector pointing out of the non-rigid material
                    int normAxis = fmobj->dflag[1];
                    if(normAxis==1 || normAxis==-1)
                    {   norm.x = (double)normAxis;
                        norm.y = 0.;
                        norm.z = 0.;
                    }
                    else if(normAxis==2 || normAxis==-2)
                    {   norm.x = 0.;
                        norm.y = normAxis>0 ? 1. : -1. ;
                        norm.z = 0.;
                    }
                    else if(normAxis==3 || normAxis==-3)
                    {   norm.x = 0.;
                        norm.y = 0.;
                        norm.z = normAxis>0 ? 1. : -1. ;
                    }
                    else
                    {   double radAngle=(double)fmobj->dflag[1]*PI_CONSTANT/180.;
                        norm.x=cos(radAngle);
                        norm.y=-sin(radAngle);
                        norm.z = 0.;
                    }
                }
                else if(fmobj->dflag[0]==4)
                {	// use special normals for cutting simulation with rake angle in dflag[1]
                    // and the material below the crack as the first defined material
                    // Assumes material 1 = material to cut (need not be used), 2 is tool, and 3 is roller bar
                    //	(note: theID is one less than the number)
                    Vector nrpos;
                    CopyScaleVector(&nrpos,&mvf[i]->disp,1./massi);
                    //if(MaterialBase::GetFieldMatID(i)==0)
                    if(nrpos.y<0. || MaterialBase::GetFieldMatID(rigidFld)==2)
                    {	norm.x=0.;
                        norm.y=1.;
                    }
                    else if(MaterialBase::GetFieldMatID(rigidFld)!=2)
                    {	double radAngle=(double)fmobj->dflag[1]*PI_CONSTANT/180.;
                        norm.x=cos(radAngle);
                        norm.y=-sin(radAngle);
                    }
                    norm.z=0.;
                }
            
                // get approach direction momentum form delPi.n (actual (vr-vi).n = delPi.n/mi)
                dotn=DotVectors(&delPi,&norm);
            
                // 3a. With this check, any movement apart will be taken as noncontact
                // Also, frictional contact assumes dotn<0
                if(dotn>=0.) break;
            
                // 3b. Displacement check
                if(contact.displacementCheck)
                {	// rigid material displacement was scaled by volume, while non-rigid was weighted by mass
					CopyScaleVector(&rigidDisp,&mvf[rigidFld]->disp,1./actualRigidVolume);
                    CopyScaleVector(&dispi,&mvf[i]->disp,1./massi);
                    hasDisplacements = TRUE;
                
                    // convert dotn to velocity of approach
                    double dvel = dotn/massi;
                
                    // check for contact
                    if(contact.MaterialContact(&dispi,&rigidDisp,&norm,dvel,postUpdate,deltime)==SEPARATED) break;
                }
                
                // passed all tests
                inContact = TRUE;
                break;
            }
            
            // continue if not in contact, unless it is an imperfect interface law
            if(!inContact && (maxContactLaw!=IMPERFECT_INTERFACE)) continue;
        }
		
		switch(maxContactLaw)
		{	case STICK:
			case NOCONTACT:
                // use delPi as is to get rigid particle velocity
                // NOCONTACT is here always
                // STICK here only if in contact
				break;
				
			case FRICTIONLESS:
				CopyScaleVector(&delPi,&norm,dotn);
				break;
				
			case FRICTIONAL:
                GetFrictionalDeltaMomentum(&delPi,&norm,dotn,maxFriction);
				break;
				
                
            case IMPERFECT_INTERFACE:
			{	// get displacement of material i and the rigid material
				if(!hasDisplacements)
				{   // rigid material displacement was scaled by volume, while non-rigid was weighted by mass
					CopyScaleVector(&rigidDisp,&mvf[rigidFld]->disp,1./actualRigidVolume);
					CopyScaleVector(&dispi,&mvf[i]->disp,1./massi);
				}
				
				// find displacement difference vector
				Vector delta=rigidDisp;
				SubVector(&delta,&dispi);
				
				// get interface forces for future use and nonuniform grid correction
				Vector fImp;
				double rawEnergy;
				
				// Get raw surface area, it is divided by hperp in GetInterfaceForces()
				// Scale voltot=voli+rigidVolume to voltot*sqrt(2*vmin/voltot)
				double rawSurfaceArea = sqrt(2.*fmin(voli,rigidVolume)*(voli+rigidVolume));
				
				// times nodal position for axisymmetric (above volumes were areas)
				if(fmobj->IsAxisymmetric()) rawSurfaceArea *= nd[nodenum]->x;
				
				bool createNode = GetInterfaceForcesForNode(&delta,&norm,Dn,Dnc,Dt,&fImp,&rawEnergy,
															rawSurfaceArea,&delPi,dotn,inContact,postUpdate,massi);
				
				if(createNode && !postUpdate)
				{	MaterialInterfaceNode::currentNode=new MaterialInterfaceNode(nd[nodenum],vfld,i,-1,&fImp,rawEnergy);
					if(MaterialInterfaceNode::currentNode==NULL)
					{	throw CommonException("Memory error allocating storage for a material interface node.",
																"CrackVelocityFieldMulti::MaterialContact");
					}
				}
                break;
			}
                
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

// Adjust change in momentum for frictional contact in tangential direction
// If has component of tangential motion, calculate force depending on whether it is sticking or sliding
void CrackVelocityFieldMulti::GetFrictionalDeltaMomentum(Vector *delPi,Vector *norm,double dotn,double frictionCoeff)
{
    // get tangential vector and its magnitude
    Vector tang;
    CopyVector(&tang,delPi);
    AddScaledVector(&tang,norm,-dotn);
    double tangMag = sqrt(DotVectors(&tang,&tang));
    
    // if has tangential motion, see if sticking or sliding
    if(!DbleEqual(tangMag,0.))
    {	ScaleVector(&tang,1./tangMag);
        double dott = DotVectors(delPi,&tang);
        
        // make it positive for comparison to the positive -frictionCoeff*dotn
        if(dott < 0.)
        {	ScaleVector(&tang,-1.);
            dott = -dott;
        }
        
        // if this is true, it is sliding and convert to frictional force
        // if not, leave delPi alone in its stick condition
        if(dott > -frictionCoeff*dotn)
        {	AddScaledVector(norm,&tang,-frictionCoeff);
            CopyScaleVector(delPi,norm,dotn);
        }
    }
}

// retrieve volume gradient
void CrackVelocityFieldMulti::GetVolumeGradient(int matfld,NodalPoint *ndptr,Vector *grad,double scale)
{
	if(fmobj->IsAxisymmetric())
	{	double nr=ndptr->x/mpmgrid.gridx;
		int n = fabs(nr)<0.01 ? 0 : 1 ;
		if(n==0)
		{	grad->x=0.;
			grad->y = mvf[matfld]->volumeGrad->y>=0. ? scale : -scale ;
		}
		else
		{	grad->x = scale*mvf[matfld]->volumeGrad->x;
			grad->y = scale*mvf[matfld]->volumeGrad->y;
		}
		grad->z=0.;
	}
	else
		CopyScaleVector(grad,mvf[matfld]->volumeGrad,scale);
}


// Handle interface forces or not (if they are perfect interfces
// Contact handled here only for perfect interface parts (Dt or Dn < 0), if done return FALSE
// Imperfect interfaces are handled next
// Use LIMIT_FORCES to decide if force is too high as determined by whether or not
//      the materials position is forced to pass the center of mass position by
//      the calculated formce
#define LIMIT_FORCES
bool CrackVelocityFieldMulti::GetInterfaceForcesForNode(Vector *delta,Vector *norm,double Dn,double Dnc,
                    double Dt,Vector *fImp,double *rawEnergy,double rawSurfaceArea,
					Vector *delPi,double dotn,bool inContact,bool postUpdate,double mred)
{
    // perpendicular distance to correct contact area and contact by positions
    double dist;
    
    // normal displacement (norm is normalized) = delta . norm, subtract adjustment when using position
	// which hae been precalculated
    double deln = DotVectors(delta,norm);
    if(!mpmgrid.GetContactByDisplacements())
	{	// Cannot use tang here, which means 3D non regular will be less accurate
		// But, all regular and all 2D will be correct without tang
		// Future Goal: get hperp without needing tang for general 3D case
		// (for efficiency, call hperp method separately)
		dist = mpmgrid.GetPerpendicularDistance(norm, NULL, 0.);
        deln -= mpmgrid.positionCutoff*dist;
	}
    
    // tangential vector in tang
    Vector tang;
    CopyVector(&tang,delta);
    AddScaledVector(&tang,norm,-deln);				// delta - deln (n) = dott (t)
    double delt=sqrt(DotVectors(&tang,&tang));
    if(!DbleEqual(delt,0.)) ScaleVector(&tang,1/delt);
	
	// if using displacements (i.e., dist not found above), find dist now with tang
	if(mpmgrid.GetContactByDisplacements())
		dist = mpmgrid.GetPerpendicularDistance(norm, &tang, delt);
    
	// initialize interfacial forces
    double trn=0.,trt=0.;
	
#ifdef LIMIT_FORCES
	// scale by minimum volume in perpendicular distance
    // Now forces are g-mm/sec^2
	double surfaceArea = rawSurfaceArea/dist;
	
    // Convert delPi to shear momentum only (delPi - dotn (n) = dott (t)), then
	//   if shear force limited leave alone, otherwise set to zero
	//   if normal force limited, add normal back, otherwise leave alone
    AddScaledVector(delPi,norm,-dotn); 
    
    // get shear (which is always be positive because sign of tang changes to accomodate it)
	// in g/(mm sec^2)
    if(Dt>=0.)
	{	// get force and compare to maximum force
		trt = 1.e6*Dt*delt*surfaceArea;
		double dott=sqrt(DotVectors(delPi,delPi));
		double maxFt = 2.*dott/timestep + 2.*mred*delt/(timestep*timestep);
        if(trt > maxFt)
		{	// limit force and retain delPi (which is for shear now)
            trt = 0.;
		}
        else
		{	// force OK, so remove shear momentum change now
            ZeroVector(delPi);
		}
    }
	// no 'else' needed because trt is zero, which means shear direction is flagged as perfect
	//     and retain delPi (which is already for shear stick)
    
	// get normal traction in g/(mm sec^2) - but different separated or in contact
	double maxFn;
	if(deln>0.)
	{	if(Dn>=0.)
		{   // normal direction in tension is imperfect
			trn = 1.e6*Dn*deln*surfaceArea;
			maxFn = 2.*dotn/timestep + 2.*mred*deln/(timestep*timestep);
			if(trn > maxFn)
			{	// limit force and add normal momentum change back into delPi
				trn = 0;
				AddScaledVector(delPi,norm,dotn);
			}
		}
		else
		{   // normal direction in tension flagged perfect
			// limit forces (trn already 0) and add normal momentum change back into delPi
			AddScaledVector(delPi,norm,dotn);
		}
	}
	else if(Dnc>=0.)
	{	// normal direction in compression is imperfect
		trn = 1.e6*Dnc*deln*surfaceArea;
		maxFn = 2.*dotn/timestep + 2.*mred*deln/(timestep*timestep);
		if(trn < maxFn)
		{	// limit negative force and add normal momentum change back into delPi
		    trn = 0;
			AddScaledVector(delPi,norm,dotn);
		}
	}
    else
    {   // normal direction is compression flagged as perfect
		// limit forces (trn already 0) and add normal momentum change back into delPi
        AddScaledVector(delPi,norm,dotn);
    }
    
    // find (trn n + trt t)*Ai for force in cartesian coordinates
    CopyScaleVector(fImp, norm, trn);
    AddScaledVector(fImp, &tang, trt);
    
	// linear elastic energy
    *rawEnergy = 0.5*(trn*deln + trt*delt);

#else
	
    if(Dt<0)
    {	if( (!inContact && Dn>=0.) || (inContact && Dnc>=0.) )
        {	// prefect in tangential, but imperfect in normal direction
            // make stick in tangential direction only by subtracting dotn.norm
            AddScaledVector(delPi,norm,-dotn);
        }
        else
        {   // else perfect in both so return with the stick conditions already in delPi
            // interface calculations will be skipped
            return FALSE;
        }
    }
	else
	{   // transverse traction in g/(mm sec^2)
        trt = 1.e6*Dt*delt;
        
		if( (!inContact && Dn<0.) || (inContact && Dnc<0.) )
		{   // perfect in normal direction, but imperfect in tangential direction
			// make stick in normal direction only
			CopyScaleVector(delPi,norm,dotn);
		}
		else
		{   // imperfect both directions, just imperfect interface forces below and nothing changed here
			// no change in momentum, just imperfect interface forces later and nothing changed here
			ZeroVector(delPi);
		}
	}
    
    if(postUpdate) return FALSE;
    
	// get normal traction in g/(mm sec^2) - but different separated or in contact
	if(deln>0.)
	{	if(Dn>=0.)
        {   // normal direction in tension is imperfect
            trn=1.e6*Dn*deln;
        }
	}
	else if(Dnc>=0.)
	{	// normal direction in compression
		trn=1.e6*Dnc*deln;
	}
    
	// scale by minimum volume in perpendicular distance
    // Now forces are g-mm/sec^2
	double surfaceArea = rawSurfaceArea/dist;
	trn *= surfaceArea;
	trt *= surfaceArea;
    
    // find (trn n + trt t)*Ai for force in cartesian coordinates
    CopyScaleVector(fImp, norm, trn);
    AddScaledVector(fImp, &tang, trt);
    
    // total energy (not increment) is (1/2)(trn dn + trt dt)*Ai in g-mm^2/sec^2
    // units will be g-mm^2/sec^2
    *rawEnergy = 0.5*(trn*deln + trt*delt);
    
#endif
	
	return TRUE;
    
}

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

// zero one component of moment and velocity
void CrackVelocityFieldMulti::SetMomVel(Vector *norm)
{	int i;
    for(i=0;i<maxMaterialFields;i++)
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

// set one component of force to -p(interpolated)/time such that updated momentum
//    of pk.i + deltime*ftot.i will be zero
void CrackVelocityFieldMulti::SetFtot(Vector *norm,double deltime)
{	int i;
    for(i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
        {	mvf[i]->SetFtotDirection(norm,deltime);
        }
    }
}

// add one component of force such that updated momentum will be mass*velocity
void CrackVelocityFieldMulti::AddFtot(Vector *norm,double deltime,double vel)
{	int i;
    for(i=0;i<maxMaterialFields;i++)
    {	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
        {	mvf[i]->AddFtotDirection(norm,deltime,vel);
        }
    }
}

#pragma mark ACCESSORS

// total number of non-rigid points
int CrackVelocityFieldMulti::GetNumberPointsNonrigid(void) { return numberPoints-numberRigidPoints; }

// total mass all velocity fields (rigid particles mass not counted)
double CrackVelocityFieldMulti::GetTotalMass(void)
{	int i;
	double mass=0;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			mass+=mvf[i]->mass;
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

// get volume for all nonrigid materials
double CrackVelocityFieldMulti::GetVolumeNonrigid(void)
{	int i;
	double volume = 0.;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveNonrigidField(mvf[i]))
			volume += mvf[i]->GetContactVolume();
	}
	return volume;
}

// get total volume for all materials
// WARNING: this doubles the volume for axisymmetric nodes at r=0
//   to enable volume screening to work
double CrackVelocityFieldMulti::GetVolumeTotal(double ndr)
{	int i;
	double volume = 0.;
	for(i=0;i<maxMaterialFields;i++)
	{	if(MatVelocityField::ActiveField(mvf[i]))
			volume += mvf[i]->GetContactVolume();
	}
	if(fmobj->IsAxisymmetric())
	{	double nr=ndr/mpmgrid.gridx;
		int n = fabs(nr)<0.01 ? 0 : 1 ;
		if(n==0) volume *= 2.;
	}
	return volume;
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
			AddVector(&fk,mvf[i]->GetFtotPtr());
	}
	return fk;
}

// add contact force on rigid material to the input vector
void CrackVelocityFieldMulti::SumAndClearRigidContactForces(Vector *fcontact,bool clearForces)
{	int rigidFieldNum;
	MatVelocityField *rigidField=GetRigidMaterialField(&rigidFieldNum);
	if(rigidField!=NULL)
	{	AddVector(fcontact,rigidField->GetFtotPtr());
		if(clearForces) ZeroVector(rigidField->GetFtotPtr());
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

/* in response to crack contact, change moment by changing velocity of all 
	nonrigid materials the same amount
 
   Change velocity by dP/M, where M is total mass
   Material i velocity becomes vi = pi/mi + dP/M
   Material i momentum change is mi vi = pi + mi dP/M
*/
void CrackVelocityFieldMulti::ChangeMomentum(Vector *delP,bool postUpdate,double deltime)
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


	
