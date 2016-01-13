/********************************************************************************
    CrackSurfaceContact.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Jan 08 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.    
********************************************************************************/

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
}

// Print contact law settings (if has one) and finalize crack law and set if has imperfect interface
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
void CrackSurfaceContact::MaterialContactPairs(int maxFields)
{
	// fill all pairs with default material properties
	mmContactLaw=(ContactLaw ***)malloc(maxFields*sizeof(ContactLaw **));
	int i,j;
	for(i=0;i<maxFields-1;i++)
	{	mmContactLaw[i]=(ContactLaw **)malloc((maxFields-1-i)*sizeof(ContactLaw *));
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

// return TRUE if imperfect interface
short CrackSurfaceContact::IsImperfect(int number) { return (short)(crackContactLaw[number]->IsImperfectInterface()); }

/*	Calculate change in momentum when there is contact. Return true or false if an adjustment was calculated
	If BC at the node, the delta momemtum should be zero in fixed direction
	Only called if both verified are verified and have 1 or more particles
	This method should ignore material that are ignoring cracks
*/
bool CrackSurfaceContact::GetDeltaMomentum(NodalPoint *np,Vector *delPa,CrackVelocityField *cva,CrackVelocityField *cvb,
											Vector *normin,int number,bool postUpdate,double deltime,int *inContact)
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
	
	// get normalized normal vector and find Delta p_a . n (actual (vb-va).n = dotn*(ma+mb)/(ma*mb))
	Vector norm;
	CopyScaleVector(&norm,normin,1./sqrt(DotVectors2D(normin,normin)));
	double dotn=DotVectors2D(delPa,&norm);
	
	// With the first check, any movement apart will be taken as noncontact
	// Also, frictional contact assume dvel<0.
	if(dotn>=0.)
		*inContact=SEPARATED;
	else
	{	// if approaching, check displacements
        // (Note: to use only velocity, skip the following displacement check)
		Vector dispa=cva->GetCMDisplacement(np,true);
		dispa.x/=massa;
		dispa.y/=massa;
		Vector dispb=cvb->GetCMDisplacement(np,true);
		dispb.x/=massb;
		dispb.y/=massb;
		
		// normal cod
		double dnorm=(dispb.x-dispa.x)*norm.x + (dispb.y-dispa.y)*norm.y
                        - mpmgrid.GetNormalCODAdjust(&norm,NULL,0);
		if(postUpdate)
		{	double dvel=(massa+massb)*dotn/(massa*massb);
			dnorm+=dvel*deltime;
		}
		
		// if current displacement positive then no contact
		if(dnorm >= 0.) *inContact=SEPARATED;
	}
	
	// if separated, then no contact unless possibly needed for an imperfect interface
	if(crackContactLaw[number]->ContactIsDone(*inContact==IN_CONTACT)) return false;
	
	// Now need to change momentum. For imperfect interface, change only for perfect directions
	double mredDE;
	if(crackContactLaw[number]->IsFrictionalContact())
	{	bool getHeating = postUpdate && ConductionTask::crackContactHeating;
		double mred = (massa*massb)/(massa+massb);
		double contactArea = 1.;
		if(crackContactLaw[number]->FrictionLawNeedsContactArea())
		{	// Angled path correction (2D only)
			double dist = mpmgrid.GetPerpendicularDistance(&norm, NULL, 0.);
			
			// Area correction method (new): sqrt(2*vmin/vtot)*vtot/dist = sqrt(2*vmin*vtot)/dist
			double vola = cva->GetVolumeNonrigid(true),volb = cvb->GetVolumeNonrigid(true),voltot=vola+volb;
			contactArea = sqrt(2.0*fmin(vola,volb)*voltot)/dist;
			if(fmobj->IsAxisymmetric()) contactArea *= np->x;
		}
		if(!crackContactLaw[number]->GetFrictionalDeltaMomentum(delPa,&norm,dotn,&mredDE,mred,
										getHeating,contactArea,*inContact==IN_CONTACT,deltime,NULL))
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
	{
		// Contact handled here only perfect interface (Dt or Dn < 0)
		// Imperfect interfaces are handled as forces later
		if(crackContactLaw[number]->IsPerfectTangentialInterface())
		{	if(!crackContactLaw[number]->IsPerfectNormalInterface(*inContact==IN_CONTACT))
			{	// prefect in tangential, but imperfect in normal direction
				// make stick in tangential direction only
				AddScaledVector(delPa,&norm,-dotn);
			}
			// else perfect in both so return with the stick conditions already in delPa
		}
		else if(crackContactLaw[number]->IsPerfectNormalInterface(*inContact==IN_CONTACT))
		{	// perfect in normal direction, but imperfect in tangential direction
			// make stick in normal direction only
			CopyScaleVector(delPa,&norm,dotn);
		}
		else
		{	// no change in momentum, just imperfect interface forces later and nothing changed here
			return false;
		}
		
	}
	
	return true;
}

// Calculate forces at imperfect interfaces and both CrackVelocityFields are present and have particles
// Return TRUE if imperfect interface or FALSE if not
// Only for cracks as imperfect interfaces
bool CrackSurfaceContact::GetInterfaceForceOnCrack(NodalPoint *np,Vector *fImp,CrackVelocityField *cva,
				CrackVelocityField *cvb,Vector *unnorm,int number,double *rawEnergy,double nodalx)
{
	// no forces needed if really perfect, was handled by contact momentum change
	if(crackContactLaw[number]->IsPerfectTangentialInterface() && crackContactLaw[number]->IsPerfectNormalInterface())
	{	return false;
	}
	
	// displacement or position
	Vector da,db;
	double mnode=1./cva->GetTotalMass(true);
	Vector dispa=cva->GetCMDisplacement(np,true);
	da.x=dispa.x*mnode;
	da.y=dispa.y*mnode;
	mnode=1./cvb->GetTotalMass(true);
	Vector dispb=cvb->GetCMDisplacement(np,true);
	db.x=dispb.x*mnode;
	db.y=dispb.y*mnode;
	
	// normal vector (assumes 2D because this is for cracks only)
	Vector norm = *unnorm;
	ScaleVector(&norm,1./sqrt(norm.x*norm.x+norm.y*norm.y));
			
    // Angled path correction (2D only)
	double dist = mpmgrid.GetPerpendicularDistance(&norm, NULL, 0.);
    
	// Area correction method (new): sqrt(2*vmin/vtot)*vtot/dist = sqrt(2*vmin*vtot)/dist
	double vola = cva->GetVolumeNonrigid(true),volb = cvb->GetVolumeNonrigid(true),voltot=vola+volb;
	double surfaceArea = sqrt(2.0*fmin(vola,volb)*voltot)/dist;
	
    // If axisymmetric, multiply by radial position (vola, volb above were areas)
    if(fmobj->IsAxisymmetric()) surfaceArea *= nodalx;
	
	// pass to imperfect interface law
	return crackContactLaw[number]->GetCrackInterfaceForce(&da,&db,&norm,surfaceArea,dist,fImp,rawEnergy);
}

// return SEPARATED if not in contact or IN_CONTACT if now in contact
// displacement is from a to b (i.e. dispbma = db-da)
// norm assummed to be normalized, dvel assumed found using normalized norm too
short CrackSurfaceContact::MaterialContact(Vector *dispbma,Vector *norm,double dvel,bool postUpdate,double deltime)
{
	// normal cod
	double dnorm=(dispbma->x*norm->x + dispbma->y*norm->y + dispbma->z*norm->z)
                    - mpmgrid.GetNormalCODAdjust(norm,NULL,0);
	
	// on post update, adjust by normal velocity difference
	if(postUpdate) dnorm+=dvel*deltime;
	
	// if current displacement positive then no contact
	return (dnorm >= 0.) ? SEPARATED : IN_CONTACT ;
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
