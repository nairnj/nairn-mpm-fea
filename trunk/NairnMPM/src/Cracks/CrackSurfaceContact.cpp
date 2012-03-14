/********************************************************************************
    CrackSurfaceContact.cpp
    NairnMPM
    
    Created by John Nairn on Thu Jan 08 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.    
********************************************************************************/

// 1. uncomment _VELOCITY_ONLY_ to use origin contact checked based only on velocity (crack contact only)
//		(don't need for multimaterial contact, instead set Vmin to 0 and Dcheck to 0)
// 2. commment out _VELOCITY_ONLY_ to check based velocity as possibly other criteria
//#define _VELOCITY_ONLY_

#include "Cracks/CrackSurfaceContact.hpp"
#include "Nodes/NodalPoint.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Materials/MaterialBase.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"

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
	Dnc=-101.;					// perfect in compression by default
	Dt=-1.;						// perfect in shear by default
	hasImperfectInterface=FALSE;	// flag for any imperfect interfaces
	moveOnlySurfaces=TRUE;			// move surfaces, plane moves at midpoint of surfaces
	preventPlaneCrosses=TRUE;		// if true, move surfaces that cross the crack plane back to the crack plane
	contactByDisplacements=TRUE;	// contact by displacements
	positionCutoff=0.8;				// element fraction when contact by positions
	materialFriction=0.;				// material contact friction
	materialDn=-1.;						// prefect in tension by default
	materialDnc=-101.;					// not set in compression by default
	materialDt=-1.;						// perfect in shear by default
	materialContactVmin=0.0;			// cutoff to kick in other contact checks
	displacementCheck=FALSE;			// if implementing check on displacement or position (last thing)
	materialNormalMethod=MAXIMUM_VOLUME_GRADIENT;		// method to find normals in multimaterial contact
	rigidGradientBias=1.;				// Use rigid gradient unless material volume gradient is this much higher (only normal method 2)
}

// Print contact law settings for cracks and finalize variables
void CrackSurfaceContact::Output(int numberOfCracks)
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
		if(Dnc<-100.) Dnc=Dn;
		sprintf(hline,"imperfect interface\n     Dnt = %g MPa/mm, Dnc = %g MPa/mm, Dt = %g MPa/mm",
				Dn,Dnc,Dt);
		hasImperfectInterface=TRUE;
	}
	else
	{   ContactLaw=FRICTIONAL;
		sprintf(hline,"frictional with coefficient of friction: %.6f",friction);
	}
	
	// print results
	cout << "Default Contact: " << hline << endl;
#ifdef _VELOCITY_ONLY_
	cout << "Contact Detection: Normal dv < 0" << endl;
#else
	cout << "Contact Detection: Normal cod < 0 AND normal dv < 0" << endl;
	if(contactByDisplacements)
		cout << "   (normal cod from displacements)" << endl;
	else
		cout << "   (normal cod from position with contact when separated less than " << positionCutoff << " of cell)" << endl;
#endif
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
		if(customDnc<-100.) customDnc=customDn;
		sprintf(hline,"imperfect interface: Dn = %g MPa/mm, Dnc = %g MPa/mm, Dt = %g MPa/mm",
				customDn,customDnc,customDt);
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
		if(materialDnc<100.) materialDnc=materialDn;
		sprintf(hline,"imperfect interface\n     Dnt = %g MPa/mm, Dnc = %g MPa/mm, Dt = %g MPa/mm",
				materialDn,materialDnc,materialDt);
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
#ifdef _VELOCITY_ONLY_
	cout << join << "(Normal dv < 0)" << endl;
#else
	cout << join << "(Normal dv < 0)";
	strcpy(join," & ");
	if(displacementCheck)
	{	cout << join << "(Normal cod < 0)" << endl;
		if(contactByDisplacements)
			cout << "   (normal cod from displacements)" << endl;
		else
			cout << "   (normal cod from position with contact when separated less than " << positionCutoff << " of cell)" << endl;
	}
	else
		cout << endl;
#endif
	
	cout << "Normal Calculation: ";
	switch(materialNormalMethod)
	{	case MAXIMUM_VOLUME_GRADIENT:
			cout << " gradient of material with maximum volume gradient," << endl;
			cout << "                        but prefer rigid material with bias factor = " << rigidGradientBias;
			rigidGradientBias*=rigidGradientBias;
			break;
		case MAXIMUM_VOLUME:
			cout << " gradient of material with maximum volume";
			break;
		case EACH_MATERIALS_MASS_GRADIENT:
			cout << " each material's own mass gradient";
			break;
		case AVERAGE_MAT_VOLUME_GRADIENTS:
			cout << " average volume gradient of material and other material with largest volume";
			break;
		default:
			break;
	}
	cout << endl;
    
    // development flags for multimaterial contact
    if(fmobj->dflag[0] > 0)
    {   cout << "** Development flag for custom contact **" << endl;
        switch(fmobj->dflag[0])
        {   case 1:
                cout << "   Each material uses its own mass gradient" << endl;
                break;
            case 2:
                cout << "   Average mass gradient of the contacting materials" << endl;
                break;
            case 3:
                cout << "   Specified normal along axis " << fmobj->dflag[1] << " with 1,2,3 for x,y,z" << endl;
            case 4:
                cout << "   Special normals for cutting. Top of tool using rake angle " << fmobj->dflag[1] <<
                ". Bottom of tool normal = (0,1)." << endl;
                break;
            default:
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
	
	// check all active materials
	for(i=0;i<nmat;i++)
	{	int mati=theMaterials[i]->GetField();
		if(mati<0) continue;
		for(j=0;j<nmat;j++)
		{	int matj=theMaterials[j]->GetField();
			if(matj<0 || i==j) continue;		// skip if no field or same material
			
			// look from custom friction from mat i to mat j
			ContactDetails *pairContact=theMaterials[i]->GetContactToMaterial(j+1);
			if(pairContact==NULL) continue;
			
			if(mati<matj)
				mmContact[mati][matj-mati-1]=*pairContact;
			else
				mmContact[matj][mati-matj-1]=*pairContact;
		}
	}
}

#pragma mark CrackSurfaceContact: Extrapolation Methods

// In task 1, track displacements or position and track volume of the entire crack velocity field
void CrackSurfaceContact::AddDisplacementVolumeTask1(short vfld,int matfld,NodalPoint *ndpt,MPMBase *mptr,double shape)
{	// exit if has no cracks and is in single material mode (i.e., not contact being done)
	if(firstCrack==NULL && maxMaterialFields==1) return;
	
	// displacement or position for contact calculations
	if(contactByDisplacements)
	{	Vector pdisp=mptr->pos;
		ndpt->AddDisplacement(vfld,matfld,mptr->mp*shape,SubVector(&pdisp,&mptr->origpos));
	}
	else
		ndpt->AddDisplacement(vfld,matfld,mptr->mp*shape,&mptr->pos);
	
	double rho=theMaterials[mptr->MatID()]->rho*0.001;	// in g/mm^3
	ndpt->AddUnscaledVolume(vfld,mptr->mp*shape/rho);
}

// In task 6, track displacements or position but do not need to track volume
void CrackSurfaceContact::AddDisplacementTask6(short vfld,int matfld,NodalPoint *ndpt,MPMBase *mptr,double shape)
{	// exit if has no cracks and is in single material mode (i.e., not contact being done)
	if(firstCrack==NULL && maxMaterialFields==1) return;
	
	// displacement or position for contact calculations
	if(contactByDisplacements)
	{	Vector pdisp=mptr->pos;
		ndpt->AddDisplacement(vfld,matfld,mptr->mp*shape,SubVector(&pdisp,&mptr->origpos));
	}
	else
		ndpt->AddDisplacement(vfld,matfld,mptr->mp*shape,&mptr->pos);
}

#pragma mark CrackSurfaceContact: Contact Calculations

// return TRUE if any contact being done
short CrackSurfaceContact::HasContact(int number) { return (short)(CrackContactLaw[number].law!=NOCONTACT); }

// return TRUE if imperfect interface
short CrackSurfaceContact::IsImperfect(int number) { return (short)(CrackContactLaw[number].law==IMPERFECT_INTERFACE); }

#ifndef _BC_CRACK_SIDE_ONLY_
/*	Calculate change in momentum when there is contact. Return TRUE or FALSE if an
	adjustment was calculated
	If BC at the node, the delta momemtum should be zero in fixed direction
	This called only if _BC_CRACK_SIDE_ONLY_ not defined in NodalPointMPM.cpp
	Only called if both verified are verified and have 1 or more particles
*/
short CrackSurfaceContact::GetDeltaMomentum(NodalPoint *np,Vector *delPa,CrackVelocityField *cva,CrackVelocityField *cvb,
											Vector *normin,int number,bool postUpdate,double deltime)
{
	// first determine if there is contact
	short inContact=IN_CONTACT;
	
	// velocities above and below
	double massa=cva->GetTotalMass();
	Vector pka=cva->GetCMatMomentum();
	double massb=cvb->GetTotalMass();
	Vector pkb=cvb->GetCMatMomentum();
	double mnode=1./(massa+massb);
	
	// screen low masses
	double aratio=massa*mnode;
	if(aratio<1.e-6 || aratio>0.999999) return FALSE;
	
	// find Delta p_a (see notes)
	CopyScaleVector(delPa,&pkb,massa*mnode);
	AddScaledVector(delPa,&pka,-massb*mnode);
	
	// get normalized normal vector and find Delta p_a . n (actual (vb-va).n = dotn*(ma+mb)/(ma*mb))
	Vector norm;
	CopyScaleVector(&norm,normin,1./sqrt(DotVectors2D(normin,normin)));
	double dotn=DotVectors2D(delPa,&norm);
	
#ifdef _VELOCITY_ONLY_
	// old version, which used to check only velocity
	if(dotn>=0.) inContact=SEPARATED;
#else
	
	// With the check, any movement apart will be taken as noncontact
	// Also, frictional contact assume dvel<0.
	if(dotn>=0.)
		inContact=SEPARATED;
	else
	{	// if approach, check displacements
		Vector dispa=cva->GetCMDisplacement();
		dispa.x/=massa;
		dispa.y/=massa;
		Vector dispb=cvb->GetCMDisplacement();
		dispb.x/=massb;
		dispb.y/=massb;
		
		// normal cod
		double dnorm=(dispb.x-dispa.x)*norm.x + (dispb.y-dispa.y)*norm.y;
		if(!contactByDisplacements) dnorm-=normalCODAdjust;
		if(postUpdate)
		{	double dvel=(massa+massb)*dotn/(massa*massb);
			dnorm+=dvel*deltime;
		}
		
		// if current displacement positive then no contact
		if(dnorm >= 0.) inContact=SEPARATED;
	}
#endif
	
	// if separated, then no contact unless possibly needed for an imperfect interface
	if(inContact==SEPARATED && CrackContactLaw[number].law!=IMPERFECT_INTERFACE) return FALSE;
	
	// Now need to change momentum. For imperfect interface, may or may not need a chnage
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
				}
			}
            break;
			
		case IMPERFECT_INTERFACE:
			// Contact handled here only perfect interface (Dt or Dn < 0)
			// Imperfect interfaces are handled as forces later
			if(CrackContactLaw[number].Dt<0)
			{	if( (inContact==SEPARATED && CrackContactLaw[number].Dn>=0.) ||
				   (inContact==IN_CONTACT && CrackContactLaw[number].Dnc>=0.) )
				{	// prefect in tangential, but imperfect in normal direction
					// make stick in tangential direction only
					AddScaledVector(delPa,&norm,-dotn);
				}
				// else perfect in both so return with the stick conditions already in delPa
			}
			else if( (inContact==SEPARATED && CrackContactLaw[number].Dn<0.) ||
					(inContact==IN_CONTACT && CrackContactLaw[number].Dnc<0.) )
			{	// perfect in normal direction, but imperfect in tangential direction
				// make stick in normal direction only
				CopyScaleVector(delPa,&norm,dotn);
			}
			else
			{	// no change in momentum, just imperfect interface forces later and nothing changed here
				return FALSE;
			}
			break;
			
        default:
            break;
    }
	
	return TRUE;
}

#else
/*	Calculate change in momentum when there is contact. Return TRUE or FALSE if an
		adjustment was calculated
	Calculate above and below crack which are equal and opposite unless there is also
		a velocity BC on the same node
	This called only if _BC_CRACK_SIDE_ONLY_ defined in NodalPointMPM.cpp
	Only called if both verified are verified and have 1 or more particles
*/
short CrackSurfaceContact::GetDeltaMomentum(NodalPoint *np,Vector *delPa,Vector *delPb,
		CrackVelocityField *cva,CrackVelocityField *cvb,Vector *norm,int number,bool postUpdate,double deltime)
{
	if this options is turned on need to check this code and include check for coontact
	also check all code within _BC_CRACK_SIDE_ONLY_ sections
	and delete CheckBCCMVelocity() which needs to be replaced by new strategy
		
	// find cm velocity and stick changes
	double massa=cva->GetTotalMass();
	double massb=cvb->GetTotalMass();
	Vector pka=cva->GetCMatMomentum();
	Vector pkb=cvb->GetCMatMomentum();
	double mnode=1./(massa+massb);
	Vector vcm;
	vcm.x=(pka.x+pkb.x)*mnode;
	vcm.y=(pka.y+pkb.y)*mnode;
	bool fixedCM=np->CheckBCCMVelocity(&vcm);
	
	// stick change above the crack: = - m(a) ( v(a) - v(cm) )
	delPa->x=massa*vcm.x-pka.x;
	delPa->y=massa*vcm.y-pka.y;
	
	// stick change below the crack: = - m(b) ( v(b) - v(cm) )
	if(fixedCM)
	{	delPb->x=massb*vcm.x-pkb.x;
		delPb->y=massb*vcm.y-pkb.y;
	}
	else
	{	delPb->x=-delPa->x;
		delPb->y=-delPa->y;
	}
	
    switch(CrackContactLaw[number].law)
    {	case STICK:
			if(inContact==SEPARATED) return FALSE;
            break;
        
        case FRICTIONLESS:
			if(inContact==SEPARATED) return FALSE;
			TangentialSlipDeltaP(delPa,norm);
			if(fixedCM)
				TangentialSlipDeltaP(delPb,norm);
			else
			{	delPb->x=-delPa->x;
				delPb->y=-delPa->y;
			}
            break;
		
        case FRICTIONAL:
			if(inContact==SEPARATED) return FALSE;
			FrictionalDeltaP(delPa,norm,number);
			if(fixedCM)
				FrictionalDeltaP(delPb,norm,number);
			else
			{	delPb->x=-delPa->x;
				delPb->y=-delPa->y;
			}
            break;
		
		case IMPERFECT_INTERFACE:
			if(CrackContactLaw[number].Dt<0)
			{	if( (inContact==SEPARATED && CrackContactLaw[number].Dn<0.) ||
						(inContact==IN_CONTACT && CrackContactLaw[number].Dnc<0.) )
				{	// both contact laws are perfect - keep as stick
					break;
				}
				else
				{	// prefect in tangential, but imperfect in normal direction
					NormalSlipDeltaP(delPa,norm);
					if(fixedCM)
						NormalSlipDeltaP(delPa,norm);
					else
					{	delPb->x=-delPa->x;
						delPb->y=-delPa->y;
					}
				}
			}
			else if( (inContact==SEPARATED && CrackContactLaw[number].Dn<0.) ||
						(inContact==IN_CONTACT && CrackContactLaw[number].Dnc<0.) )
			{	// perfect in normal direction, but imperfect in tangential direction
				TangentialSlipDeltaP(delPa,norm);
				if(fixedCM)
					TangentialSlipDeltaP(delPb,norm);
				else
				{	delPb->x=-delPa->x;
					delPb->y=-delPa->y;
				}
			}
			else
			{	// no change in momentum, just imperfect interface forces later
				return FALSE;
			}
			break;
			
		default:
			break;
	}
	
	return TRUE;
}
#endif

// find frictionaless tangnential slip change in momentum where on input
//   delP is stick change in momentum
//   norm is unnormalized normal into material (but sign is irrelevant)
void CrackSurfaceContact::TangentialSlipDeltaP(Vector *delP,Vector *norm)
{
	// Find mass times changes in normal direction velocity
	// mdelvn = - m ( v - v(cm) ) . n = delP . n (now normalized)
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
	// mdelvt = - m ( v - v(cm) ) . t = delP . t (now normalized)
    double mdelvt=(delP->x*norm->y - delP->y*norm->x)/(norm->x*norm->x + norm->y*norm->y);
	delP->x=mdelvt*norm->y;
	delP->y=-mdelvt*norm->x;
}

// find frictional change in momentum where on input
//   delP is stick change in momentum
//   norm is unnormalized normal into material (but sign is irrelevant)
void CrackSurfaceContact::FrictionalDeltaP(Vector *delP,Vector *norm,int number)
{
	// mdelvn = - m ( v - v(cm) ) . n = delP . n (unnormalized)
	double mdelvn=delP->x*norm->x + delP->y*norm->y;
	
	// mdelvt = - m ( v - v(cm) ) . t = delP . t (unnormalized)
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
short CrackSurfaceContact::GetInterfaceForce(NodalPoint *np,Vector *fImp,CrackVelocityField *cva,
				CrackVelocityField *cvb,Vector *norm,int number,double *rawEnergy)
{
	// no forces needed if really perfect, was handled by contact momentum change
	if(CrackContactLaw[number].Dn<0. && CrackContactLaw[number].Dt<0. && CrackContactLaw[number].Dnc<0.)
		return FALSE;
	
	// displacement or position
	Vector da,db;
	double mnode=1./cva->GetTotalMass();
	Vector dispa=cva->GetCMDisplacement();
	da.x=dispa.x*mnode;
	da.y=dispa.y*mnode;
	mnode=1./cvb->GetTotalMass();
	Vector dispb=cvb->GetCMDisplacement();
	db.x=dispb.x*mnode;
	db.y=dispb.y*mnode;
			
	double dnunnorm,dtunnorm,trn=0.,trt=0.;
	
	if(CrackContactLaw[number].Dn>=0. || CrackContactLaw[number].Dnc>=0.)
	{	// normal displacement (unnormalized - missing 1/|norm|)
		dnunnorm=(db.x-da.x)*norm->x + (db.y-da.y)*norm->y;
		if(!contactByDisplacements) dnunnorm-=sqrt(norm->x*norm->x+norm->y*norm->y)*normalCODAdjust;
		// Normal traction in g/(mm sec^2) - but different separated or in contact
		if(dnunnorm>0.)
		{	// normal direction in tension
			if(CrackContactLaw[number].Dn>=0.)
				trn=1.e6*CrackContactLaw[number].Dn*dnunnorm;
			else
			{	// interface perfect in tension, if also perfect in shear can exit
				dnunnorm=0.;
				if(CrackContactLaw[number].Dt<0.) return FALSE;
			}
		}
		else
		{	// normal direction in compression
			if(CrackContactLaw[number].Dnc>=0.)
				trn=1.e6*CrackContactLaw[number].Dnc*dnunnorm;
			else
			{	// interface perfect in compression, if also perfect in shear can exit
				dnunnorm=0.;
				if(CrackContactLaw[number].Dt<0.) return FALSE;
			}
		}
	}
	else
	{	// perfect in normal direction
		dnunnorm=0.;
	}
			
	if(CrackContactLaw[number].Dt>=0.)
	{	// transverse force (unnormalized - missing 1/|norm|)
		dtunnorm=(db.x-da.x)*norm->y - (db.y-da.y)*norm->x;
		// transverse traction in g/(mm sec^2)
		trt=1.e6*CrackContactLaw[number].Dt*dtunnorm;
	}
	else
	{	// perfect in normal direction
		dtunnorm=0.;
	}
			
	// find trn n + trt t and finally normalize
	double norm2=1./(norm->x*norm->x + norm->y*norm->y);	// square because |norm| for trn or trt and |norm| for norm->x,y
	fImp->x=(trn*norm->x+trt*norm->y)*norm2;
	fImp->y=(trn*norm->y-trt*norm->x)*norm2;
	
	// total energy (not increment) is (1/2)(trn dnunnorm + trt dtunnorm)/(norm2*norm2) in g/sec^2
	// Use norm2 because |norm| for trn and trt and |norm| for dnunnorm and dtunnorm
	// units wiill be g/sec^2
	*rawEnergy=(trn*dnunnorm + trt*dtunnorm)*norm2/2.;
	
	return TRUE;
}

// return SEPARATED if not in contact or IN_CONTACT if now in contact
// displacement is from a to b (i.e. db-da)
// norm assummed to be normalized, dvel assumed found using noramlized norm too
short CrackSurfaceContact::MaterialContact(Vector *dispa,Vector *dispb,Vector *norm,double dvel,bool postUpdate,double deltime)
{
	// normal cod
	double dnorm=((dispb->x-dispa->x)*norm->x + (dispb->y-dispa->y)*norm->y + (dispb->z-dispa->z)*norm->z);
	if(!contactByDisplacements) dnorm-=normalCODAdjust;
	
	// on post update, adjust by normal velocity difference
	if(postUpdate) dnorm+=dvel*deltime;
	
	// if current displacement positive then no contact
	return (dnorm >= 0.) ? SEPARATED : IN_CONTACT ;
}

#pragma mark ACCESSORS

// return current setting for moving only surfaces
bool CrackSurfaceContact::GetMoveOnlySurfaces(void) { return moveOnlySurfaces; }
void CrackSurfaceContact::SetMoveOnlySurfaces(bool moveSurfaces) { moveOnlySurfaces=moveSurfaces; }

// return current setting for moving only surfaces
bool CrackSurfaceContact::GetPreventPlaneCrosses(void) { return preventPlaneCrosses; }
void CrackSurfaceContact::SetPreventPlaneCrosses(bool preventCross) { preventPlaneCrosses=preventCross; }

// return current setting for contact method
bool CrackSurfaceContact::GetContactByDisplacements(void) { return contactByDisplacements; }
void CrackSurfaceContact::SetContactByDisplacements(bool newContact) { contactByDisplacements=newContact; }

// convert fraction to actual position offset for contact by positions
void CrackSurfaceContact::SetNormalCODCutoff(double meshSize) { normalCODAdjust=positionCutoff*meshSize; }
double CrackSurfaceContact::GetNormalCODCutoff(void) { return normalCODAdjust; }

// material contact law for field mati to field matj
short CrackSurfaceContact::GetMaterialContactLaw(int mati,int matj)
{	// index based on smaller of the two indices
	return mati<matj ? mmContact[mati][matj-mati-1].law : mmContact[matj][mati-matj-1].law ;
}

// material coefficient of friction for field mati to field matj
double CrackSurfaceContact::GetMaterialFriction(int mati,int matj)
{	// index based on smaller of the two indices
	return mati<matj ? mmContact[mati][matj-mati-1].friction : mmContact[matj][mati-matj-1].friction ;
}
