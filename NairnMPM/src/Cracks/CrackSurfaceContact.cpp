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
#include "Exceptions/MPMWarnings.hpp"
#include "Nodes/CrackVelocityFieldMulti.hpp"
#include "Nodes/MaterialContactNode.hpp"

// Single global contact law object
CrackSurfaceContact crackContact;

extern double timestep;

#pragma mark CrackSurfaceContact: Constructors and Destructors

// Constructors
CrackSurfaceContact::CrackSurfaceContact()
{
	moveOnlySurfaces = true;			// move surfaces, plane moves at midpoint of surfaces
	preventPlaneCrosses = false;		// if true, move surfaces that cross the crack plane back to the crack plane
	
	// contact
	contactLawID = -1;
	contactByDisplacements = true;		// contact by displacements for cracks
}

// Print contact law settings for cracks and finalize variables
// throws std::bad_alloc, CommonException()
void CrackSurfaceContact::Output(void)
{
	// allocate memory for custom crack contact laws
	char *p=new char[(numberOfCracks+1)*sizeof(ContactLaw *)];
	contactLaw =( ContactLaw **)p;
	
	// Global crack contact law (must be set, if not force to frictionless)
	contactLawID = MaterialBase::GetContactLawNum(contactLawID);
	if(contactLawID<0)
		throw CommonException("Crack settings must select a default contact law","CrackSurfaceContact::Output");
	contactLaw[0] = (ContactLaw *)theMaterials[contactLawID];
	cout << "Default Contact Law: " << contactLaw[0]->name << " (number " << (contactLawID+1) << ")" << endl;
	if(contactLaw[0]->IsImperfectInterface()) mpmgrid.hasImperfectInterface=true;
	
	// print other settings
	cout << "Contact Detection: Normal cod < 0 AND normal dv < 0" << endl;
    mpmgrid.OutputContactByDisplacements(false,contactByDisplacements,positionCutoff);
    
    cout << "Crack Particle Position Updates: ";
    cout << "XPIC(1)" << endl;
	if(GetMoveOnlySurfaces())
		cout << "Crack Plane Updating: Average of the crack surfaces" << endl;
	else
		cout << "Crack Plane Updating: Use center of mass velocity" << endl;
	if(GetPreventPlaneCrosses())
		cout << "Crack Plane Crosses: surface particles moved back to the current plane" << endl;
	else
		cout << "Crack Plane Crosses: ignored" << endl;
	
	cout << "Crack Shape Functions: Classic" << endl;
}

// Print contact law settings (if has one) and finalize crack law and see if has imperfect interface
// throws CommonException()
void CrackSurfaceContact::CustomCrackContactOutput(int &customCrackID,int number)
{
	// no custom law was set, so pick the default law
	if(customCrackID<0)
	{	contactLaw[number] = contactLaw[0];
		return;
	}
	
	// custom law
	customCrackID = MaterialBase::GetContactLawNum(customCrackID);
	if(customCrackID<0)
		throw CommonException("Custom crack contact must select a valid contact law","CrackSurfaceContact::Output");
	contactLaw[number] = (ContactLaw *)theMaterials[customCrackID];
	cout << "    Custom Contact Law: " << contactLaw[number]->name << " (number " << (customCrackID+1) << ")" << endl;
	if(contactLaw[number]->IsImperfectInterface()) mpmgrid.hasImperfectInterface = true;
}

#pragma mark CrackSurfaceContact: Contact Calculations

// return TRUE if any contact being done
short CrackSurfaceContact::HasContact(int number) { return (short)(!contactLaw[number]->IgnoreContact()); }

/*	Calculate change in momentum when there is contact. Return true or false if an adjustment was calculated
	If BC at the node, the delta momemtum should be zero in fixed direction
	Only called if both velocities are verified and have 1 or more particles
	This method should ignore materials that are ignoring cracks
*/
bool CrackSurfaceContact::GetDeltaMomentum(NodalPoint *np,Vector *delPa,CrackVelocityField *cva,CrackVelocityField *cvb,
								Vector *normin,int number,int callType,double deltime,int *inContact,FMPMContact *ccCache)
{
#pragma mark Crack Contact: Get values and Delta Pa
	// velocities above and below
	bool hasParticles;
	double massa,massb;
	Vector pka=cva->GetCMatMomentum(hasParticles,&massa,NULL,false);
	Vector pkb=cvb->GetCMatMomentum(hasParticles,&massb,NULL,false);
	double mnode=1./(massa+massb);
	
	// screen low masses
	if(contactLaw[number]->IsImperfectInterface())
	{	double aratio=massa*mnode;
		if(aratio<1.e-6 || aratio>0.999999)
		{	MaterialContactNode::ContactSetLowMass(ccCache,true,1);
			return false;
		}
	}
	
	// find Delta p_a (see notes)
	CopyScaleVector(delPa,&pkb,massa*mnode);
	AddScaledVector(delPa,&pka,-massb*mnode);
	
#pragma mark Crack Contact: Get normal, normal separation deltaDotn and vector deltaDotn
	// get normalized normal vector and find dPDotn = Delta p_a . n (actual (vb-va).n = dPDotn*(ma+mb)/(ma*mb))
	Vector norm;
	ZeroVector(&norm); // for 2D
	CopyScaleVector(&norm,normin,1./sqrt(DotVectors(normin,normin)));
	double dPDotn = DotVectors(delPa,&norm);
	
	// get displacements (dispa, diapb) and cod (delta)
	Vector dispa,dispb;
	if(contactLaw[number]->IsImperfectInterface())
	{	// imperfect interface overrides contactByDicplacement to always use
		// displacements; i.e. interface should start in contact
		dispa = cva->GetCMDisplacement(np,true,true);
		dispb = cvb->GetCMDisplacement(np,true,true);
	}
	else
	{	// contact law displacement calculations - uses contactByDisplacement if set
		dispa = cva->GetCMDisplacement(np,true,contactByDisplacements);
		dispb = cvb->GetCMDisplacement(np,true,contactByDisplacements);
	}
	ScaleVector(&dispa, (1. / massa));
	ScaleVector(&dispb, (1. / massb));
	Vector delta = dispb;
	SubVector(&delta,&dispa);
	
	// Get normal displacement (interfaces don't use this deltaDotn, so contactByDisplacements doesn't matter)
	double deltaDotn = MaterialSeparation(DotVectors(&dispb,&norm),DotVectors(&dispa,&norm),&norm,np,
								   contactByDisplacements,positionCutoff);
	
	// on post update, adjust by normal velocity difference
	double mred = massa*massb*mnode;
	if(callType == UPDATE_MOMENTUM_CALL)
	{	double dvel = dPDotn/mred;
		deltaDotn += dvel*deltime;
	}

	// contact area
	double contactArea = 1.;
	if(contactLaw[number]->ContactLawNeedsContactArea())
	{	double vola = cva->GetContactVolumeNonrigid(true);
		double volb = cvb->GetContactVolumeNonrigid(true);
		contactArea = CrackVelocityFieldMulti::GetContactArea(np,vola,volb,&norm,NULL);
	}

	// first determine if there is contact
	*inContact=IN_CONTACT;

#pragma mark Crack Contact: Implement contact law
	// Paramters calculated above
	// np: node be analyzed
	// number: crack number
	// delPa: momentum change side a to com velocity
	// delta: separation vector (only needed by interfaces)
	// norm: normal vector
	// deltaDotn: delta.n or normal separation between materials
	// dPDotn: delPa.n
	// masses: massa, massb, mnode=1/(massa+massb)
	// contactLaw[number]: contact law
	// contactArea: area only if needed
	// mredDelWf, getHeating

	// cache contact info for FMPM
	MaterialContactNode::SetContactInfo(ccCache,deltaDotn,1,delPa,&norm,contactLaw[number],contactArea);

	// Now need to change momentum. For imperfect interface, change only for perfect directions
	double mredDelWf;
	if(contactLaw[number]->IsFrictionalContact())
	{	// Handle frictional contact
		
		// second order heating needs acceleration term
		// Find (ma Fb - mb Fa)/(ma+mb)
		bool getHeating = false;
		Vector delFi;
		Vector *delFiPtr = NULL;
		if(callType==UPDATE_MOMENTUM_CALL)
		{	Vector Fb = cvb->GetCMatFtot();
			Vector Fa = cva->GetCMatFtot();
			CopyScaleVector(&delFi,&Fb,massa*mnode);
			AddScaledVector(&delFi,&Fa,-massb*mnode);
			delFiPtr = &delFi;
			if(ccCache!=NULL) ccCache->startDelFi = delFi;

			// second order heating needs acceleration too
			if(ConductionTask::matContactHeating) getHeating = true;
		}
		
		if(!contactLaw[number]->GetFrictionalDeltaMomentum(delPa,&norm,dPDotn,deltaDotn,&mredDelWf,mred,
										getHeating,contactArea,deltime,delFiPtr,np,true))
		{	*inContact = SEPARATED;
			return false;
		}
		if(mredDelWf>0.)
		{	double qrate = mredDelWf/mred;
			// As heat source need nJ/sec or multiply by 1/timestep
			// Note that this is after transport rates are calculated (by true in last parameter)
			conduction->AddFluxCondition(np,fabs(qrate/deltime),true);
			
			// this line seems to crack XCode analyzer, comment our temporarily to analyze code
#pragma omp atomic
			NodalPoint::frictionWork += qrate;
			
			if(ccCache!=NULL) ccCache->priorMredDelWf = mredDelWf;
		}
	}
	else
	{	// Handle imperfect interfaces
		
		// get tangDel (unit tangent vecotr) and deltaDott (= delta.t)
		Vector tangDel;
		double deltaDott = CrackVelocityFieldMulti::GetTangentCOD(&norm,&delta,&tangDel,ccCache);
		
		// Find delFi = (ma Fb - mb Fa)/Mc (when needed)
		bool postUpdate = callType==UPDATE_MOMENTUM_CALL;
		Vector fImp;
		if(postUpdate)
		{	Vector Fb = cvb->GetCMatFtot();
			Vector Fa = cva->GetCMatFtot();
			CopyScaleVector(&fImp,&Fb,massa*mnode);
			AddScaledVector(&fImp,&Fa,-massb*mnode);
			if(ccCache!=NULL) ccCache->startDelFi = fImp;
		}
		else
			ZeroVector(&fImp);
		
		// Get interface force, energy, and possible changed delPa
		double rawEnergy;
		contactLaw[number]->GetInterfaceForces(&norm,&fImp,&rawEnergy,contactArea,delPa,
												dPDotn,mred,&tangDel,deltaDotn,deltaDott,postUpdate);
		
		if(postUpdate)
		{	// add force (if any) to momentum change
			AddScaledVector(delPa, &fImp, timestep);
			
			// Add interface energy. (Legacy units g-mm^2/sec^2 or multiply by 1e-9 to get J - kg-m^2/sec^2)
#pragma omp atomic
			NodalPoint::interfaceEnergy += rawEnergy;
		}
		
		// If no interface force, then should stick with returned momentum
		// If has force, still may need to change using returned altered momentum (although it could be zero)
	}
	
	return true;
}

// Return separation in normal direction
// If extraplated position, correct for grid effects by two methods
// dispb.n, dispa.n, norm, and ndptr (norm assummed to be normalized)
double CrackSurfaceContact::MaterialSeparation(double dbdotn,double dadotn,Vector *norm,NodalPoint *ndptr,bool useDisps,double cutoff)
{
	// get dnorm and correct if needed
	double dnorm;
	if(useDisps)
		dnorm = dbdotn-dadotn;
	else
	{	double r = cutoff;
		Vector dist = mpmgrid.GetPerpendicularDistance(norm,ndptr);
		if(r>0.)
		{	dnorm = dbdotn-dadotn - r*dist.x;
		}
		else
		{	// power law fit from nanoindentation paper to improve linear cutoff method
			r = -r;
			double pa = dadotn - ndptr->x*norm->x - ndptr->y*norm->y - ndptr->z*norm->z;
			double da = pa>0. ? 2.*pow(pa/(1.25*dist.x),r) - 1. : 1 - 2.*pow(-pa/(1.25*dist.x),r);
			double pb = dbdotn - ndptr->x*norm->x - ndptr->y*norm->y - ndptr->z*norm->z;
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

