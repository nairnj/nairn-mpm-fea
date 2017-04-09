/********************************************************************************
    CrackSurfaceContact.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Jan 08 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.    
********************************************************************************/

#include "stdafx.h"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Nodes/NodalPoint.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Materials/MaterialBase.hpp"
#include "Materials/ContactLaw.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "System/UnitsController.hpp"
#include "Exceptions/CommonException.hpp"

// Single global contact law object
CrackSurfaceContact contact;

extern double timestep;

#pragma mark CrackSurfaceContact: Constructors and Destructors

// Constructors
CrackSurfaceContact::CrackSurfaceContact()
{
	crackContactLawID=-1;
	
	hasImperfectInterface=false;	// flag for any imperfect interfaces
	moveOnlySurfaces=true;			// move surfaces, plane moves at midpoint of surfaces
	preventPlaneCrosses=false;		// if true, move surfaces that cross the crack plane back to the crack plane
	
	materialContactLawID=-1;
	
	materialContactVmin=0.0;			// cutoff to kick in other contact checks
	displacementCheck=true;			// if implementing check on displacement or position (last thing)
	materialNormalMethod=AVERAGE_MAT_VOLUME_GRADIENTS;		// method to find normals in multimaterial contact
	rigidGradientBias=1.;				// Use rigid gradient unless material volume gradient is this much higher (only normal method 2)
}

// Print contact law settings for cracks and finalize variables
// throws std::bad_alloc, CommonException()
void CrackSurfaceContact::Output(void)
{
	// allocate memory for custom crack contact laws
	char *p=new char[(numberOfCracks+1)*sizeof(ContactLaw *)];
	crackContactLaw=(ContactLaw **)p;
	
	// Global material contact law (must be set,if not force to frictionless)
	crackContactLawID = MaterialBase::GetContactLawNum(crackContactLawID);
	if(crackContactLawID<0)
		throw CommonException("Crack settings must select a default contact law","CrackSurfaceContact::Output");
	crackContactLaw[0] = (ContactLaw *)theMaterials[crackContactLawID];
	cout << "Default Contact Law: " << crackContactLaw[0]->name << " (number " << (crackContactLawID+1) << ")" << endl;
	if(crackContactLaw[0]->IsImperfectInterface()) hasImperfectInterface=true;
	
	// print other settings
	cout << "Contact Detection: Normal cod < 0 AND normal dv < 0" << endl;
    mpmgrid.OutputContactByDisplacements();
	if(GetMoveOnlySurfaces())
		cout << "Crack Plane Updating: Average of the crack surfaces" << endl;
	else
		cout << "Crack Plane Updating: Use center of mass velocity" << endl;
	if(GetPreventPlaneCrosses())
		cout << "Crack Plane Crosses: surface particles moved back to the current plane" << endl;
	else
		cout << "Crack Plane Crosses: ignored" << endl;
	
#ifdef CRACK_GIMP
	if(mpmgrid.crackParticleSize<0.)
		mpmgrid.crackParticleSize = 0.;
	else if(mpmgrid.crackParticleSize>1.)
		mpmgrid.crackParticleSize = 1.;
	cout << "Crack particle size for shape functions: " << mpmgrid.crackParticleSize << endl;
#endif
}

// Print contact law settings (if has one) and finalize crack law and set if has imperfect interface
// throws CommonException()
void CrackSurfaceContact::CustomCrackContactOutput(int &customCrackID,int number)
{
	// no custom law was set
	if(customCrackID<0)
	{	crackContactLaw[number] = crackContactLaw[0];
		return;
	}
	
	// custom law
	customCrackID = MaterialBase::GetContactLawNum(customCrackID);
	if(customCrackID<0)
		throw CommonException("Custom crack contact must select a default contact law","CrackSurfaceContact::Output");
	crackContactLaw[number] = (ContactLaw *)theMaterials[customCrackID];
	cout << "    Custom Contact Law: " << crackContactLaw[number]->name << " (number " << (customCrackID+1) << ")" << endl;
	if(crackContactLaw[number]->IsImperfectInterface()) hasImperfectInterface=true;
}

// Print contact law settings for cracks and finalize variables
// throws CommonException()
void CrackSurfaceContact::MaterialOutput(void)
{
	// Global material contact law (must be set,if not force to frictionless)
	materialContactLawID = MaterialBase::GetContactLawNum(materialContactLawID);
	if(materialContactLawID<0)
		throw CommonException("Multimaterial mode must select a default contact law","CrackSurfaceContact::MaterialOutput");
	materialContactLaw = (ContactLaw *)theMaterials[materialContactLawID];
	cout << "Default Contact Law: " << materialContactLaw->name << " (number " << (materialContactLawID+1) << ")" << endl;
	if(materialContactLaw->IsImperfectInterface()) hasImperfectInterface=true;
	
	// print contact detection method
	char join[3];
	join[0]=0;
	cout << "Contact Detection: ";
	if(materialContactVmin>0.)
	{	cout << "(Vrel >= " << materialContactVmin << ")";
		strcpy(join," & ");
	}
	cout << join << "(Normal dv < 0)";
	strcpy(join," & ");
	if(displacementCheck)
	{	cout << join << "(Normal cod < 0)" << endl;
        mpmgrid.OutputContactByDisplacements();
	}
	else
		cout << endl;
	
	cout << "Normal Calculation: ";
	switch(materialNormalMethod)
	{	case MAXIMUM_VOLUME_GRADIENT:
			cout << " gradient of material or paired material (if all nonrigid), or the rigid material (if" << endl;
			cout << "                   one rigid material), that has highest magnitude. When has rigid" << endl;
			cout << "                   material, prefer rigid material gradient with bias factor = " << rigidGradientBias;
			rigidGradientBias*=rigidGradientBias;       // squared as needed in algorithm
			break;
		case MAXIMUM_VOLUME:
			cout << " gradient of material with maximum volume";
			break;
		case AVERAGE_MAT_VOLUME_GRADIENTS:
			cout << " volume-weighted mean gradient of material and other materials lumped (if all nonrigid)," << endl;
			cout << "                   on just the rigid material (if one rigid material). When has rigid" << endl;
			cout << "                   material, prefer rigid material gradient with bias factor = " << rigidGradientBias;
			rigidGradientBias*=rigidGradientBias;       // squared as needed in algorithm
			break;
		case EACH_MATERIALS_MASS_GRADIENT:
			cout << " each material's own mass gradient";
			break;
		case SPECIFIED_NORMAL:
			cout << " use the specified normal of ";
			PrintVector("",&contactNormal);
		default:
			break;
	}
	cout << endl;
}

// prepare array for material contact details
// throws CommonException()
void CrackSurfaceContact::MaterialContactPairs(int maxFields)
{
	// create double array of pairs
	mmContactLaw = new (nothrow) ContactLaw **[maxFields];
	if(mmContactLaw==NULL)
	{	throw CommonException("Memory error creating contact pairs array","CrackSurfaceContact::MaterialContactPairs");
	}
	
	// fill all pairs with default material properties
	int i,j;
	for(i=0;i<maxFields-1;i++)
	{	mmContactLaw[i] = new (nothrow) ContactLaw *[maxFields-1-i];
		if(mmContactLaw[i]==NULL)
		{	throw CommonException("Memory error creating contact pairs array","CrackSurfaceContact::MaterialContactPairs");
		}
		
		// to default law
		for(j=i+1;j<maxFields;j++) mmContactLaw[i][j-i-1] = materialContactLaw;
	}
	
	// check all active materials and change laws that were specified
	for(i=0;i<nmat;i++)
	{	int mati=theMaterials[i]->GetField();			// may be a shared field
		if(mati<0) continue;							// skip if not used
		
		// loop over all other materials
		for(j=0;j<nmat;j++)
		{	int matj=theMaterials[j]->GetField();		// may be a shared field
			if(matj<0 || i==j) continue;				// skip if no field or same material
			
			// look from custom friction from mat i to mat j
			int pairContactID=theMaterials[i]->GetContactToMaterial(j+1);
			if(pairContactID<0) continue;
			pairContactID = MaterialBase::GetContactLawNum(pairContactID);
			
			// setting more than one shared material overwrite previous ones
			if(mati<matj)
			{	mmContactLaw[mati][matj-mati-1]=(ContactLaw *)theMaterials[pairContactID];
				if(mmContactLaw[mati][matj-mati-1]->IsImperfectInterface()) hasImperfectInterface=true;
			}
			else
			{	mmContactLaw[matj][mati-matj-1]=(ContactLaw *)theMaterials[pairContactID];
				if(mmContactLaw[matj][mati-matj-1]->IsImperfectInterface()) hasImperfectInterface=true;
			}
		}
	}
}

#pragma mark CrackSurfaceContact: Contact Calculations

// return TRUE if any contact being done
short CrackSurfaceContact::HasContact(int number) { return (short)(!crackContactLaw[number]->IgnoreContact()); }

/*	Calculate change in momentum when there is contact. Return true or false if an adjustment was calculated
	If BC at the node, the delta momemtum should be zero in fixed direction
	Only called if both verified are verified and have 1 or more particles
	This method should ignore material that are ignoring cracks
*/
bool CrackSurfaceContact::GetDeltaMomentum(NodalPoint *np,Vector *delPa,CrackVelocityField *cva,CrackVelocityField *cvb,
											Vector *normin,int number,int callType,double deltime,int *inContact)
{
	// first determine if there is contact
	*inContact=IN_CONTACT;
	
	// velocities above and below
	bool hasParticles;
	double massa,massb;
	Vector pka=cva->GetCMatMomentum(hasParticles,&massa);
	Vector pkb=cvb->GetCMatMomentum(hasParticles,&massb);
	double mnode=1./(massa+massb);
	
	// screen low masses
	double aratio=massa*mnode;
	if(aratio<1.e-6 || aratio>0.999999) return false;
	//if(aratio<1.e-3 || aratio>0.999) return FALSE;
	
	// find Delta p_a (see notes)
	CopyScaleVector(delPa,&pkb,massa*mnode);
	AddScaledVector(delPa,&pka,-massb*mnode);
	
	// get normalized normal vector and find dPDotn = Delta p_a . n (actual (vb-va).n = dPDotn*(ma+mb)/(ma*mb))
	Vector norm;
	CopyScaleVector(&norm,normin,1./sqrt(DotVectors2D(normin,normin)));
	double dPDotn = DotVectors2D(delPa,&norm);
	
	// will need to get displacements if doing displacement check or
	// if contact is imperfect interface
	double deltaDotn=0.;
	Vector delta = MakeVector(0.,0.,0.);
	if(dPDotn<0. || crackContactLaw[number]->IsImperfectInterface())
	{	// displacement calculations
		Vector dispa=cva->GetCMDisplacement(np,true);
		dispa.x/=massa;
		dispa.y/=massa;
		dispa.z = 0.;
		Vector dispb=cvb->GetCMDisplacement(np,true);
		dispb.x/=massb;
		dispb.y/=massb;
		dispb.z = 0.;
		delta = dispb;
		SubVector(&delta,&dispa);
		deltaDotn = MaterialSeparation(&delta,&dispa,&dispb,&norm,np);
	}
	
	// With the first check, any movement apart will be taken as noncontact
	// Also, frictional contact assume dvel<0.
	if(dPDotn >= 0.)
		*inContact=SEPARATED;
	else
	{	// on post update, adjust by normal velocity difference
		if(callType!=MASS_MOMENTUM_CALL)
		{	double dvel=(massa+massb)*dPDotn/(massa*massb);
			deltaDotn += dvel*deltime;
		}
		
		// if current displacement positive then no contact
		*inContact = (deltaDotn >= 0.) ? SEPARATED : IN_CONTACT ;
	}
	
	// if separated, then no contact unless possibly needed for an imperfect interface
	if(crackContactLaw[number]->ContactIsDone(*inContact==IN_CONTACT)) return false;
	
	// Now need to change momentum. For imperfect interface, change only for perfect directions
	double mredDE;
	double mred = (massa*massb)/(massa+massb);
	if(crackContactLaw[number]->IsFrictionalContact())
	{	bool getHeating = (callType==UPDATE_MOMENTUM_CALL) && ConductionTask::crackContactHeating;
		double contactArea = 1.;
		if(crackContactLaw[number]->ContactLawNeedsContactArea())
		{	// Angled path correction (cracks are only 2D)
			Vector dist = mpmgrid.GetPerpendicularDistance(&norm,np);
			
			// Area correction method (new): sqrt(2*vmin/vtot)*vtot/dist = sqrt(2*vmin*vtot)/dist
			// dist weightings to allow for Tartan grid
			double vola = cva->GetVolumeNonrigid(true),volb = cvb->GetVolumeNonrigid(true),voltot=vola+volb;
			contactArea = sqrt(2.0*fmin(vola*dist.y,volb*dist.z)*voltot)/dist.x;
			if(fmobj->IsAxisymmetric()) contactArea *= np->x;
		}
		
		// second order heating needs acceleration term
		// Find (ma Fb - mb Fa)/(ma+mb)
		Vector at;
		Vector *atPtr = NULL;
		if(getHeating)
		{	Vector Fb = cvb->GetCMatFtot();
			Vector Fa = cva->GetCMatFtot();
			CopyScaleVector(&at,&Fb,massa*mnode);
			AddScaledVector(&at,&Fa,-massb*mnode);
			atPtr = &at;
		}
		
		if(!crackContactLaw[number]->GetFrictionalDeltaMomentum(delPa,&norm,dPDotn,&mredDE,mred,
										getHeating,contactArea,*inContact==IN_CONTACT,deltime,atPtr))
		{	return false;
		}
		if(mredDE>0.)
		{	double qrate = mredDE/mred;
			NodalPoint::frictionWork += qrate;
		
			// As heat source need nJ/sec or multiply by 1/timestep
			// Note that this is after transport rates are calculated (by true in last parameter)
			conduction->AddFluxCondition(np,fabs(qrate/deltime),true);
		}
	}
	else
	{	// get tangDel and deltaDott
		Vector tangDel;
		double deln = DotVectors(&delta,&norm);			// delta.n, but not same as correct deltaDotn from above
		CopyVector(&tangDel,&delta);
		AddScaledVector(&tangDel,&norm,-deln);				// delta - deln (n) = deltaDott (t)
		
		// by normalizing to positive delt, hat t always points in positive direction
		double deltaDott=sqrt(DotVectors(&tangDel,&tangDel));
		if(!DbleEqual(deltaDott,0.)) ScaleVector(&tangDel,1./deltaDott);
		
		// get contact area - angled path correction (cracks are only 2D)
		Vector dist = mpmgrid.GetPerpendicularDistance(&norm,np);
		
		// Area correction method (new): sqrt(2*vmin/vtot)*vtot/dist = sqrt(2*vmin*vtot)/dist
		// dist weightings to allow for Tartan grid
		double vola = cva->GetVolumeNonrigid(true),volb = cvb->GetVolumeNonrigid(true),voltot=vola+volb;
		double contactArea = sqrt(2.0*fmin(vola*dist.y,volb*dist.z)*voltot)/dist.x;
		if(fmobj->IsAxisymmetric()) contactArea *= np->x;
		
		// get input force if needed and then call for interface force and energy
		// Find delFi = (ma Fb - mb Fa)/Mc (when needed)
		Vector fImp;
		if(callType==UPDATE_MOMENTUM_CALL)
		{	Vector Fb = cvb->GetCMatFtot();
			Vector Fa = cva->GetCMatFtot();
			CopyScaleVector(&fImp,&Fb,massa*mnode);
			AddScaledVector(&fImp,&Fa,-massb*mnode);
		}
		else
			ZeroVector(&fImp);
		double rawEnergy;
		crackContactLaw[number]->GetInterfaceForces(&norm,&fImp,&rawEnergy,
													contactArea,delPa,dPDotn,mred,&tangDel,deltaDotn,deltaDott,dist.x);
		
#ifndef MANDMIMPINT
		if(callType==UPDATE_MOMENTUM_CALL)
		{	// add force (if any) to momentum change
			AddScaledVector(delPa, &fImp, timestep);
			
			// Add interface energy. (Legacy units g-mm^2/sec^2 or multiply by 1e-9 to get J - kg-m^2/sec^2)
			NodalPoint::interfaceEnergy += rawEnergy;
		}
#else
		if(callType==MASS_MOMENTUM_CALL)
		{	// add only prior to update and add to forces
			cva->AddFtotSpreadTask3(&fImp);
			ScaleVector(&fImp,-1.);
			cvb->AddFtotSpreadTask3(&fImp);
			
			// Add interface energy. (Legacy units g-mm^2/sec^2 or multiply by 1e-9 to get J - kg-m^2/sec^2)
			NodalPoint::interfaceEnergy += rawEnergy;
		}
#endif
		
		// If no interface force, then should stick with returned momentum
		// If has force, still may need to change using returned altered momentum (although it could be zero)
	}
	
	return true;
}

// Return dispbma.n
// If extraplated position correct for edge effects
// Input displacement from material a and material b, their difference, norm, and ndptr
//		norm assummed to be normalized
double CrackSurfaceContact::MaterialSeparation(Vector *dispbma,Vector *dispa,Vector *dispb,Vector *norm,NodalPoint *ndptr)
{
	// get dnorm and correct if needed
	double dnorm;
	if(mpmgrid.GetContactByDisplacements())
		dnorm = dispbma->x*norm->x + dispbma->y*norm->y + dispbma->z*norm->z;
	else
	{	double r = mpmgrid.positionCutoff;
		Vector dist = mpmgrid.GetPerpendicularDistance(norm,ndptr);
		if(r>0.)
		{	dnorm = dispbma->x*norm->x + dispbma->y*norm->y + dispbma->z*norm->z;
			dnorm -= r * dist.x;
		}
		else
		{	r = -r;
			double pa = (dispa->x-ndptr->x)*norm->x + (dispa->y-ndptr->y)*norm->y + (dispa->z-ndptr->z)*norm->z;
			double da = pa>0. ? 2.*pow(pa/(1.25*dist.x),r) - 1. : 1 - 2.*pow(-pa/(1.25*dist.x),r);
			double pb = (dispb->x-ndptr->x)*norm->x + (dispb->y-ndptr->y)*norm->y + (dispb->z-ndptr->z)*norm->z;
			double db = pb>0. ? 2.*pow(pb/(1.25*dist.x),r) - 1. : 1 - 2.*pow(-pb/(1.25*dist.x),r);
			dnorm = (db-da)*dist.x;
			//cout << "# (" << ndptr->x << "," << ndptr->y << ") " << pa << "," << pb << "," << da << "," << db << "," << dnorm << endl;
		}
	}
	
	return dnorm;
}

#pragma mark ACCESSORS

// return current setting for moving only surfaces
bool CrackSurfaceContact::GetMoveOnlySurfaces(void) const { return moveOnlySurfaces; }
void CrackSurfaceContact::SetMoveOnlySurfaces(bool moveSurfaces) { moveOnlySurfaces=moveSurfaces; }

// return current setting for moving only surfaces
bool CrackSurfaceContact::GetPreventPlaneCrosses(void) const { return preventPlaneCrosses; }
void CrackSurfaceContact::SetPreventPlaneCrosses(bool preventCross) { preventPlaneCrosses=preventCross; }

// material contact law for field mati to field matj
ContactLaw *CrackSurfaceContact::GetMaterialContactLaw(int mati,int matj)
{	// index based on smaller of the two indices
	return mati<matj ? mmContactLaw[mati][matj-mati-1] : mmContactLaw[matj][mati-matj-1] ;
}

// set contact normal when normal is specified
void CrackSurfaceContact::SetContactNormal(double polarAngle,double aximuthAngle)
{
	double angle,sinp;
	
	if(fmobj->IsThreeD())
	{	angle = PI_CONSTANT*polarAngle/180.;
		contactNormal.z = cos(angle);
		sinp = sin(angle);
	}
	else
	{	contactNormal.z = 0.;
		sinp = 1.0;
	}
	angle = PI_CONSTANT*aximuthAngle/180.;
	contactNormal.x = cos(angle)*sinp;
	contactNormal.y = sin(angle)*sinp;
}
