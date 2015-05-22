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
#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "System/UnitsController.hpp"

// Single global contact law object
CrackSurfaceContact contact;

extern double timestep;

#pragma mark CrackSurfaceContact: Constructors and Destructors

// Constructors
CrackSurfaceContact::CrackSurfaceContact()
{
	ContactLaw=FRICTIONLESS;	// the law
	friction=0.;				// crack contact friction
	Dn=-1.;						// prefect in tension by default
	Dnc=-101.e6;					// <-100e6 means not set and should be set same as Dn
	Dt=-1.;						// perfect in shear by default
	hasImperfectInterface=FALSE;	// flag for any imperfect interfaces
	moveOnlySurfaces=TRUE;			// move surfaces, plane moves at midpoint of surfaces
	preventPlaneCrosses=FALSE;		// if true, move surfaces that cross the crack plane back to the crack plane
	materialFriction=0.;				// material contact friction
	materialDn=-1.;						// prefect in tension by default
	materialDnc=-101.e6;					// <-100e6 means not set and should be set same as Dn
	materialDt=-1.;						// perfect in shear by default
	materialContactVmin=0.0;			// cutoff to kick in other contact checks
	displacementCheck=TRUE;			// if implementing check on displacement or position (last thing)
	materialNormalMethod=AVERAGE_MAT_VOLUME_GRADIENTS;		// method to find normals in multimaterial contact
	rigidGradientBias=1.;				// Use rigid gradient unless material volume gradient is this much higher (only normal method 2)
}

// Print contact law settings for cracks and finalize variables
void CrackSurfaceContact::Output(void)
{
	char hline[200];
	
	// Default contact law
	if(friction<-10.)
	{   ContactLaw=NOCONTACT;
		sprintf(hline,"contacts ignored");
	}
	else if(friction<-.5)
	{   ContactLaw=STICK;
		sprintf(hline,"stick conditions");
	}
	else if(DbleEqual(friction,(double)0.))
	{   ContactLaw=FRICTIONLESS;
		sprintf(hline,"frictionless sliding");
	}
	else if(friction>10.)
	{   ContactLaw=IMPERFECT_INTERFACE;
		if(Dnc<-100.e6) Dnc=Dn;
		const char *label = UnitsController::Label(INTERFACEPARAM_UNITS);
		sprintf(hline,"imperfect interface\n     Dnt = %g %s, Dnc = %g %s, Dt = %g %s",
				Dn*UnitsController::Scaling(1.e-6),label,
				Dnc*UnitsController::Scaling(1.e-6),label,
				Dt*UnitsController::Scaling(1.e-6),label);
		hasImperfectInterface=TRUE;
	}
	else
	{   ContactLaw=FRICTIONAL;
		sprintf(hline,"frictional with coefficient of friction: %.6f",friction);
	}
	
	// print results
	cout << "Default Contact: " << hline << endl;
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
    
	// allocate memory for custom crack contact laws
	char *p=new char[(numberOfCracks+1)*sizeof(ContactDetails)];
    
    // this is the default contact law
	CrackContactLaw=(ContactDetails *)p;
	CrackContactLaw[0].law=ContactLaw;
	CrackContactLaw[0].friction=friction;
	CrackContactLaw[0].Dn=Dn;
	CrackContactLaw[0].Dnc=Dnc;
	CrackContactLaw[0].Dt=Dt;
}

// Print contact law settings and finalize variables
void CrackSurfaceContact::CrackOutput(bool custom,double customFriction,double customDn,double customDnc,
									  double customDt,int number)
{
	if(!custom || (customFriction<10. && DbleEqual(friction,customFriction)) || 
	   (friction>10. && customFriction>10. && DbleEqual(Dn,customDn) && DbleEqual(Dt,customDt) &&  DbleEqual(Dnc,customDnc)))
	{	CrackContactLaw[number].law=ContactLaw;
		CrackContactLaw[number].friction=friction;
		CrackContactLaw[number].Dn=Dn;
		CrackContactLaw[number].Dnc=Dnc;
		CrackContactLaw[number].Dt=Dt;
		return;
	}
	
	char hline[200];
	
	CrackContactLaw[number].friction=customFriction;
	CrackContactLaw[number].Dn=customDn;
	CrackContactLaw[number].Dnc=customDnc;
	CrackContactLaw[number].Dt=customDt;
	if(customFriction<-10.)
	{   CrackContactLaw[number].law=NOCONTACT;
		sprintf(hline,"contacts ignored");
	}
	else if(customFriction<-.5)
	{   CrackContactLaw[number].law=STICK;
		sprintf(hline,"stick conditions");
	}
	else if(customFriction>10.)
	{   CrackContactLaw[number].law=IMPERFECT_INTERFACE;
		if(customDnc<-100.e6) customDnc=customDn;
		const char *label = UnitsController::Label(INTERFACEPARAM_UNITS);
		sprintf(hline,"imperfect interface: Dn = %g %s, Dnc = %g %s, Dt = %g %s",
				customDn*UnitsController::Scaling(1.e-6),label,
				customDnc*UnitsController::Scaling(1.e-6),label,
				customDt*UnitsController::Scaling(1.e-6),label);
		hasImperfectInterface=TRUE;
	}
	else if(customFriction>0.)
	{   CrackContactLaw[number].law=FRICTIONAL;
		sprintf(hline,"frictional with coefficient of friction: %.6f",customFriction);
	}
	else
	{   CrackContactLaw[number].law=FRICTIONLESS;
		CrackContactLaw[number].friction=0.;		// to be sure
		sprintf(hline,"frictionless sliding");
	}
	cout << "    Custom Contact: " << hline << endl;
}

// Print contact law settings for cracks and finalize variables
void CrackSurfaceContact::MaterialOutput(void)
{
	char hline[200];
	
	// Global material contact
	if(materialFriction<-10.)
	{   materialContactLaw=NOCONTACT;
		sprintf(hline,"contact nodes revert to center of mass velocity field");
	}
	else if(materialFriction<-.5)
	{   materialContactLaw=STICK;
		sprintf(hline,"stick conditions");
	}
	else if(DbleEqual(materialFriction,(double)0.))
	{   materialContactLaw=FRICTIONLESS;
		sprintf(hline,"frictionless sliding");
	}
	else if(materialFriction>10.)
	{   materialContactLaw=IMPERFECT_INTERFACE;
		if(materialDnc<-100.e6) materialDnc=materialDn;
		const char *label = UnitsController::Label(INTERFACEPARAM_UNITS);
		sprintf(hline,"imperfect interface\n     Dnt = %g %s, Dnc = %g %s, Dt = %g %s",
				materialDn*UnitsController::Scaling(1.e-6),label,
				materialDnc*UnitsController::Scaling(1.e-6),label,
				materialDt*UnitsController::Scaling(1.e-6),label);
	}
	else
	{   materialContactLaw=FRICTIONAL;
		sprintf(hline,"frictional with coefficient of friction: %.6f",materialFriction);
	}
	
	// print results
	char join[3];
	join[0]=0;
	cout << "Default Contact: " << hline << endl;
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
			cout << endl;
		default:
			break;
	}
	cout << endl;
    
    // development flags for multimaterial contact
    if(fmobj->dflag[0] > 0)
    {   cout << "** Development flag for custom contact **" << endl;
        switch(fmobj->dflag[0])
		{	case 4:
                cout << "   Special normals for cutting. Top of tool using ";
                if(fmobj->dflag[1]>-90.)
                    cout << "rake angle " << fmobj->dflag[1];
                else
                    cout << "calculated normals";
                cout << ". Bottom of tool normal = (0,1)." << endl;
                break;
			case 5:
				cout << "   Radial normal for spherical inclusion" <<endl;
				break;
            default:
				cout << "   Unknown, or no longer implemented, custom contact option" << endl;
                break;
        }
    }
	
}

// prepare array for material contact details
void CrackSurfaceContact::MaterialContactPairs(int maxFields)
{
	// fill all pairs with default material properties
	mmContact=(ContactDetails **)malloc(maxFields*sizeof(ContactDetails *));
	int i,j;
	for(i=0;i<maxFields-1;i++)
	{	mmContact[i]=(ContactDetails *)malloc((maxFields-1-i)*sizeof(ContactDetails));
		for(j=i+1;j<maxFields;j++)
		{	mmContact[i][j-i-1].law=materialContactLaw;
			mmContact[i][j-i-1].friction=materialFriction;
			mmContact[i][j-i-1].Dn=materialDn;
			mmContact[i][j-i-1].Dnc=materialDnc;
			mmContact[i][j-i-1].Dt=materialDt;
		}
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
			ContactDetails *pairContact=theMaterials[i]->GetContactToMaterial(j+1);
			if(pairContact==NULL) continue;
			
			// setting more than one shared material overwrite previous ones
			if(mati<matj)
				mmContact[mati][matj-mati-1]=*pairContact;
			else
				mmContact[matj][mati-matj-1]=*pairContact;
		}
	}
}

#pragma mark CrackSurfaceContact: Contact Calculations

// return TRUE if any contact being done
short CrackSurfaceContact::HasContact(int number) { return (short)(CrackContactLaw[number].law!=NOCONTACT); }

// return TRUE if imperfect interface
short CrackSurfaceContact::IsImperfect(int number) { return (short)(CrackContactLaw[number].law==IMPERFECT_INTERFACE); }

/*	Calculate change in momentum when there is contact. Return true or false if an adjustment was calculated
	If BC at the node, the delta momemtum should be zero in fixed direction
	Only called if both verified are verified and have 1 or more particles
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
	if(*inContact==SEPARATED && CrackContactLaw[number].law!=IMPERFECT_INTERFACE) return false;
	
	// Now need to change momentum. For imperfect interface, may or may not need a change
	Vector tang;
	double dott,mu;
	
    switch(CrackContactLaw[number].law)
    {	case STICK:
            break;
			
        case FRICTIONLESS:
			CopyScaleVector(delPa,&norm,dotn);
            break;
			
        case FRICTIONAL:
			CopyVector(&tang,delPa);
			AddScaledVector(&tang,&norm,-dotn);
			dott=sqrt(DotVectors2D(&tang,&tang));
			if(!DbleEqual(dott,0.))
			{	ScaleVector(&tang,1./dott);
				dott=DotVectors2D(delPa,&tang);
				if(dott<0.)
				{	ScaleVector(&tang,-1.);
					dott=-dott;
				}
				mu=-CrackContactLaw[number].friction;
				if(dott>mu*dotn)
				{	AddScaledVector(&norm,&tang,mu);
					CopyScaleVector(delPa,&norm,dotn);
                    
                    // get frictional heating part - this is g mm^2/sec^2 = nJ
                    // Note: only add frictional heating during momentum update (when frictional
                    //   force is appropriate) and only if conduction is on.
                    if(postUpdate && conduction && ConductionTask::crackContactHeating)
                    {   if(np->NodeHasNonrigidParticles())
                        {   Vector Ftdt;
                            CopyScaleVector(&Ftdt,&tang,mu*dotn);
                            double qrate = (massa+massb)*DotVectors2D(&Ftdt,delPa)/(massa*massb);
                            
                            // As heat source need nJ/sec or multiply by 1/timestep
                            // Note that this is after transport rates are calculated (by true in last parameter)
                            conduction->AddFluxCondition(np,fabs(qrate/deltime),true);
                        }
                    }
				}
			}
            break;
			
		case IMPERFECT_INTERFACE:
			// Contact handled here only perfect interface (Dt or Dn < 0)
			// Imperfect interfaces are handled as forces later
			if(CrackContactLaw[number].Dt<0)
			{	if( (*inContact==SEPARATED && CrackContactLaw[number].Dn>=0.) ||
				   (*inContact==IN_CONTACT && CrackContactLaw[number].Dnc>=0.) )
				{	// prefect in tangential, but imperfect in normal direction
					// make stick in tangential direction only
					AddScaledVector(delPa,&norm,-dotn);
				}
				// else perfect in both so return with the stick conditions already in delPa
			}
			else if( (*inContact==SEPARATED && CrackContactLaw[number].Dn<0.) ||
					(*inContact==IN_CONTACT && CrackContactLaw[number].Dnc<0.) )
			{	// perfect in normal direction, but imperfect in tangential direction
				// make stick in normal direction only
				CopyScaleVector(delPa,&norm,dotn);
			}
			else
			{	// no change in momentum, just imperfect interface forces later and nothing changed here
				return false;
			}
			break;
			
        default:
            break;
    }
	
	return true;
}

// find frictionaless tangnential slip change in momentum where on input
//   delP is stick change in momentum
//   norm is unnormalized normal into material (but sign is irrelevant)
void CrackSurfaceContact::TangentialSlipDeltaP(Vector *delP,Vector *norm)
{
	// Find mass times changes in normal direction velocity
	// mdelvn = - m ( v - v(ctr mass) ) . n = delP . n (now normalized)
    double mdelvn=(delP->x*norm->x + delP->y*norm->y)/(norm->x*norm->x + norm->y*norm->y);
	delP->x=mdelvn*norm->x;
	delP->y=mdelvn*norm->y;
}
// find frictionless normal slip change in momentum where on input
//   delP is stick change in momentum
//   norm is unnormalized normal into material (but sign is irrelevant)
void CrackSurfaceContact::NormalSlipDeltaP(Vector *delP,Vector *norm)
{
	// Find mass times changes in normal direction velocity
	// mdelvt = - m ( v - v(ctr mass) ) . t = delP . t (now normalized)
    double mdelvt=(delP->x*norm->y - delP->y*norm->x)/(norm->x*norm->x + norm->y*norm->y);
	delP->x=mdelvt*norm->y;
	delP->y=-mdelvt*norm->x;
}

// find frictional change in momentum where on input
//   delP is stick change in momentum
//   norm is unnormalized normal into material (but sign is irrelevant)
void CrackSurfaceContact::FrictionalDeltaP(Vector *delP,Vector *norm,int number)
{
	// mdelvn = - m ( v - v(ctr mass) ) . n = delP . n (unnormalized)
	double mdelvn=delP->x*norm->x + delP->y*norm->y;
	
	// mdelvt = - m ( v - v(ctr mass) ) . t = delP . t (unnormalized)
    double mdelvt=delP->x*norm->y - delP->y*norm->x;
	
	// relative tangentical to normal stick forces is the ratio
	double ratio=mdelvt/mdelvn;
	double mu;
	if(ratio<0.)
	{	if(CrackContactLaw[number].friction<-ratio)
			mu=-CrackContactLaw[number].friction;
		else
			return;			// return to use stick conditions
	}
	else
	{	if(CrackContactLaw[number].friction<ratio)
			mu=CrackContactLaw[number].friction;
		else
			return;			// return to use stick conditions
	}
	
	// normalize and get the change
	mdelvn/=(norm->x*norm->x + norm->y*norm->y);
	delP->x=mdelvn*(norm->x + mu*norm->y);
	delP->y=mdelvn*(norm->y - mu*norm->x);
}
	
// Calculate forces at imperfect interfaces and both CrackVelocityFields are present and have particles
// Return TRUE if imperfect interface or FALSE if not
// Only for cracks as imperfect interfaces
short CrackSurfaceContact::GetInterfaceForceOnCrack(NodalPoint *np,Vector *fImp,CrackVelocityField *cva,
				CrackVelocityField *cvb,Vector *unnorm,int number,double *rawEnergy,double nodalx)
{
	// no forces needed if really perfect, was handled by contact momentum change
	if(CrackContactLaw[number].Dn<0. && CrackContactLaw[number].Dt<0. && CrackContactLaw[number].Dnc<0.)
		return FALSE;
	
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
			
    // Angled path correction (special case because norm is not normalized (2D only)
	double dist = mpmgrid.GetPerpendicularDistance(&norm, NULL, 0.);
    
	// Area correction method (new): sqrt(2*vmin/vtot)*vtot/dist = sqrt(2*vmin*vtot)/dist
	double vola = cva->GetVolumeNonrigid(true),volb = cvb->GetVolumeNonrigid(true),voltot=vola+volb;
	double surfaceArea = sqrt(2.0*fmin(vola,volb)*voltot)/dist;
	
    // If axisymmetric, multiply by radial position (vola, volb above were areas)
    if(fmobj->IsAxisymmetric()) surfaceArea *= nodalx;
	
	double dn,dt,trn = 0.,trt = 0.;
	
	if(CrackContactLaw[number].Dn>=0. || CrackContactLaw[number].Dnc>=0.)
	{	// normal displacement
		dn = (db.x-da.x)*norm.x + (db.y-da.y)*norm.y;
		if(!mpmgrid.GetContactByDisplacements())
		{	// for efficiency used calculated dist
			dn -= mpmgrid.positionCutoff*dist;
		}
		
		// Normal traction in g/(mm sec^2) - but different separated or in contact
		if(dn>0.)
		{	// normal direction in tension
			if(CrackContactLaw[number].Dn>=0.)
				trn = CrackContactLaw[number].Dn*dn*surfaceArea;
			else
			{	// interface perfect in tension, if also perfect in shear can exit
				if(CrackContactLaw[number].Dt<0.) return FALSE;
				dn = 0.;
			}
		}
		else
		{	// normal direction in compression
			if(CrackContactLaw[number].Dnc>=0.)
				trn = CrackContactLaw[number].Dnc*dn*surfaceArea;
			else
			{	// interface perfect in compression, if also perfect in shear can exit
				if(CrackContactLaw[number].Dt<0.) return FALSE;
				dn = 0.;
			}
		}
	}
	else
	{	// perfect in normal direction
		dn = 0.;
	}
			
	if(CrackContactLaw[number].Dt>=0.)
	{	// transverse force
		dt = (db.x-da.x)*norm.y - (db.y-da.y)*norm.x;
		
		// transverse traction in g/(mm sec^2)
		trt = CrackContactLaw[number].Dt*dt*surfaceArea;
	}
	else
	{	// perfect in normal direction
		dt = 0.;
	}
			
	// find trn n + trt t and finally normalize
	fImp->x = trn*norm.x + trt*norm.y;
	fImp->y = trn*norm.y - trt*norm.x;
	
	// total energy (not increment) is (1/2)(trn dnunnorm + trt dtunnorm)/(norm2*norm2) in g/sec^2
	// Use norm2 because |norm| for trn and trt and |norm| for dnunnorm and dtunnorm
	// units wiill be g/sec^2
	*rawEnergy = (trn*dn + trt*dt)/2.;
	
	return TRUE;
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
int CrackSurfaceContact::GetMaterialContactLaw(int mati,int matj)
{	// index based on smaller of the two indices
	return mati<matj ? mmContact[mati][matj-mati-1].law : mmContact[matj][mati-matj-1].law ;
}

// material coefficient of friction for field mati to field matj
double CrackSurfaceContact::GetMaterialFriction(int mati,int matj)
{	// index based on smaller of the two indices
	return mati<matj ? mmContact[mati][matj-mati-1].friction : mmContact[matj][mati-matj-1].friction ;
}

void CrackSurfaceContact::GetMaterialInterface(int mati,int matj,double *Dn,double *Dnc,double *Dt)
{   // index based on smaller of the two indices
    int i,j;
    if(mati < matj)
    {   i = mati;
        j = matj-mati-1;
    }
    else
    {   i = matj;
        j = mati-matj-1;
    }
    *Dn = mmContact[i][j].Dn;
    *Dnc = mmContact[i][j].Dnc;
    *Dt = mmContact[i][j].Dt;
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
