/********************************************************************************
    CrackSegment.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Apr 3 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Cracks/CrackSegment.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/TractionLaw.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "System/ArchiveData.hpp"
#include "Cracks/CrackHeader.hpp"
#include "System/UnitsController.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"

extern char *app;

#pragma mark CrackSegment: Constructors and Destructors

// Constructors
CrackSegment::CrackSegment()
{
}

CrackSegment::CrackSegment(Vector *xp,int tip,int matid)
{
	cp=*xp;
	surf[0]=cp;
	surf[1]=cp;
	orig=cp;
    planeInElem=surfInElem[0]=surfInElem[1]=0;
    nextSeg=NULL;
	prevSeg=NULL;
    parent=NULL;
    
	// need to track velocities for PIC updates
	ZeroVector(&cpVel);
	ZeroVector(&surfVel[0]);
	ZeroVector(&surfVel[1]);
	
	ZeroVector(&Jint);
	ZeroVector(&sif);
	ZeroVector(&tract);
    ZeroVector(&czmdG);
    czmdG.z=-1.;
    tipMatnum=tip;
	propagationJ=0.;
    steadyState=STATIONARY;
	heating=FALSE;
	SetMatID(matid);
	historyData=NULL;
}

// when add particle, find plane element and set surface element to the same one
void CrackSegment::FindInitialElement(void)
{
	planeInElem = FindElement();
	surfInElem[0] = planeInElem;
	surfInElem[1] = planeInElem;
}

#pragma mark CrackSegment: Methods

// find current element (1 based) or return 0 if no element
// Only needed to initialize and move cracks. At only times, get
//   element ID from planeElemID()
int CrackSegment::FindElement(void)
{
    // check current element
    if(planeInElem>0)
    {	if(theElements[planeElemID()]->PtInElement(cp))
            return planeInElem;
    }
	
	// use grid code
	try
	{	planeInElem = mpmgrid.FindElementFromPoint(&cp,NULL);
	}
	catch (...)
	{	planeInElem = 0;
	}
	
    return planeInElem;
}

// find current element (1 based) or return 0 if no element for crack surface
// side = ABOVE_CRACK (1) or BELOW)_CRACK (2)
// Only needed to initialize and move cracks. At only times, get
//   element ID from surfElemID(side)
int CrackSegment::FindElement(short side)
{
    int j=side-1;
	Vector cpt;
	cpt = surf[j];

    // check current element
    if(surfInElem[j]>0)
    {	if(theElements[surfaceElemID(side)]->PtInElement(cpt))
            return surfInElem[j];
    }
	
	// use grid code
	try
	{	surfInElem[j] = mpmgrid.FindElementFromPoint(&cpt,NULL);
	}
	catch (...)
	{	surfInElem[j] = 0;
	}
	
	return surfInElem[j];
}

// Reset crack plane position from surfaces (2D) (in mm)
// only used when contact.GetMoveOnlySurfaces() is true and thus get crack position
//		from previous movement of the surfaces
void CrackSegment::MovePosition(void)
{	cp.x += dPlane.x;
	cp.y += dPlane.y;
}

// Reset crack plane position from surfaces (2D) (in mm)
void CrackSegment::MovePositionToMidpoint(void)
{	cp.x = 0.5*(surf[0].x+surf[1].x);
	cp.y = 0.5*(surf[0].y+surf[1].y);
}

// Move crack plane position (2D) (in mm)
// only used when contact.GetMoveOnlySurfaces() is FALSE and thus need to move crack plane
void CrackSegment::MovePosition(Vector *velnp1,Vector *cpAcc,double dt,double shapeNorm)
{
	// XPIC(1) position change is dX = dt*(1.5*svelnp1-0.5*prevVel-cpAcc*dt)/shapeNorm
	Vector delX = SetScaledVector(&cpVel,-0.5);
	AddScaledVector(&delX,velnp1,1.5);
	AddScaledVector(&delX,cpAcc,-dt);
    ScaleVector(&delX,dt/shapeNorm);
	
	// add to current postion
	AddVector(&cp,&delX);
	
	// update velocity to new velocity
	CopyScaleVector(&cpVel,velnp1,1./shapeNorm);
}

// Move a surface position (2D) (in mm) - must move ABOVE_CRACK and then BELOW_CRACK
// side = ABOVE_CRACK (1) or BELOW)_CRACK (2)
// Cracks use PIC update such that new velocity is svelnp1 but
// 		position change is dt*(1.5*svelnp1-surfAccc*dt-0.5*prevVel)
bool CrackSegment::MoveSurfacePosition(short side,Vector *svelnp1,Vector *surfAcc,double dt,bool hasNodes)
{
	short j=side-1;
	bool movedOther=false;
	
	// if hasNodes new velocity is in svelnp1, get increment
	Vector delX;
	if(hasNodes)
	{
        // position change is dX = dt*(1.5*svelnp1-0.5*prevVel-surfAcc*dt)
		delX = SetScaledVector(&surfVel[j],-0.5);
		AddScaledVector(&delX,svelnp1,1.5);
		AddScaledVector(&delX,surfAcc,-dt);
		ScaleVector(&delX,dt);
	}
	
	if(side==ABOVE_CRACK)
	{	// above crack is first
		if(hasNodes)
		{	dPlane = delX;					// save until below is done next
			AddVector(&surf[j],&delX);
			surfVel[j] = *svelnp1;
		}
		else
		{	ZeroVector(&dPlane);
			ZeroVector(&surfVel[j]);
		}
		hadAboveNodes = hasNodes;
	}
	else
	{	// below crack is second
		// ... if has nodes than get average movement or just this one if top did not move
		//			also move top along if it had no nodes
		// ... else if no nodes, move bottom along with the top field (if it had one)
		if(hasNodes)
		{	AddVector(&surf[j],&delX);
			surfVel[j] = *svelnp1;
			
			if(hadAboveNodes)
			{	// average movement of the two surfaces
				dPlane = SetSumVectors(&dPlane,&delX);
				ScaleVector(&dPlane,0.5);
			}
			else
			{	dPlane = delX;
				AddVector(&surf[ABOVE_CRACK-1],&dPlane);	// move above by delX from below
				surfVel[ABOVE_CRACK-1] = SetScaledVector(&delX,1./dt);
                
                // in case above moved elements now (this tells to check on return)
				movedOther = true;
			}
		}
		else if(hadAboveNodes)
		{   // only had nodes above the crack, move bottom by delX from above
            AddVector(&surf[j],&dPlane);
			surfVel[j] = SetScaledVector(&dPlane,1./dt);
		}
		else
		{	// neither side has nodes
			ZeroVector(&surfVel[j]);
		}
	}
	
	return movedOther;
}

// calculate tractions due to this segment
void CrackSegment::AddTractionForceSeg(CrackHeader *theCrack)
{	// exit if no traction law
	if(MatID()<0) return;

	// Add force on this segement to both sides of the crack
	AddTractionForceSegSide(theCrack,ABOVE_CRACK,1.);
	AddTractionForceSegSide(theCrack,BELOW_CRACK,-1.);
}

// calculate tractions on one side of crack for this segment
// add forces to material velocity fields on one side of the crack
// side = ABOVE_CRACK (1) or BELOW_CRACK (2)
void CrackSegment::AddTractionForceSegSide(CrackHeader *theCrack,int side,double scale)
{
#ifdef CONST_ARRAYS
	int nds[MAX_SHAPE_NODES];
	double fn[MAX_SHAPE_NODES];
	short cvfld[MAX_SHAPE_NODES];
#else
	int nds[maxShapeNodes];
	double fn[maxShapeNodes];
	short cvfld[maxShapeNodes];
#endif
	double fnorm = 0.;
	NodalPoint *ndi;
	int	cnum=theCrack->GetNumber();
	
	// get element and shape function to extrapolate to the node
	int js = side-1;
	int iel = surfaceElemID(side);
	theElements[iel]->GetShapeFunctionsForCracks(fn,nds,surf[js]);
	int numnds = nds[0];

#define MASS_WEIGHTED_FORCE
    bool hasNodes = false;
	// find mass from nodes seen by this crack particle
	for(int i=1;i<=numnds;i++)
	{	// Get velocity field to use
		ndi = nd[nds[i]];
		cvfld[i] = ndi->GetFieldForSurfaceParticle(side,cnum,this,false);
		
		// if has particles (and they see cracks), get weighted shape function
		if (cvfld[i] >= 0)
		{
#ifdef MASS_WEIGHTED_FORCE
            // The theory: here fnorm = sum(mi*Sip) is interpreted as mass on the crack surface particle
            // The crack surface acceleration is T/fnorm where T is stored traction
            // Extrapolate this acceleration to node ai = Sip*T/fnorm
            // Convert to force on node fi = mi*Sip*T/fnorm = (fn[i]/fnorm)*T
            fn[i] *= ndi->cvf[cvfld[i]]->GetTotalMass(true);
			fnorm += fn[i];
#else
            // The theory here is node fi = Sip*T/fnorm where fnorm = sum(Sip)
            fnorm += fn[i];
#endif
            hasNodes = true;
		}
	}
	
    // skip if none
    if(!hasNodes) return;
    
	// Get 1/fnorm once for speed and add scale to handle crack side
	fnorm = scale/fnorm;
	
	// loop over all nodes seen by this crack surface particle
	for(int i=1;i<=numnds;i++)
	{	// if has particles (and they see cracks), add force
		if(cvfld[i]>=0)
		{	ndi = nd[nds[i]];
			ndi->AddFtotSpreadTask3(cvfld[i],FTract(fnorm*fn[i]));
		}
	}
}

/* This method is called for use with traction force and energy and does several tasks
	1. Get normal pointing from below to above
 	2. Get tangential that is perpendicular to normal and in direction of tangential cod
 			delta = (delta.n)n + (delta.t)t
 			tangential is in direction of delta - (delta.n)n
 	3. If theCrack!=NULL find area for traction of this crack point (returned value)
 	4. Return nCod = delta.n and nCod = delta.t
 	5. Might be better to use splines
 In 2D
 	a. Find tangential = (tx,ty,0) from line segments arround point
 	b. Normal = (-ty,tx,0)
 in 3D
    This method overridden in CrackPoint
*/
double CrackSegment::GetNormalAndTangent(CrackHeader *theCrack,Vector *norm,Vector *tang,double &nCod,double &tCod) const
{
	// area to return
	double area = 1.;
	
	// get cod
	Vector cod = SetDiffVectors(&surf[0],&surf[1]);

	// 2D Calculations
	if(prevSeg==NULL)
	{	// changed to half in revision 1850
		tang->x=(nextSeg->cp.x-cp.x)/2.;
		tang->y=(nextSeg->cp.y-cp.y)/2.;
	}
	else if(nextSeg==NULL)
	{	// changed to half in revision 1850
		tang->x=(cp.x-prevSeg->cp.x)/2.;
		tang->y=(cp.y-prevSeg->cp.y)/2.;
	}
	else
	{	tang->x=(nextSeg->cp.x-prevSeg->cp.x)/2.;
		tang->y=(nextSeg->cp.y-prevSeg->cp.y)/2.;
	}
	tang->z = 0.;
	
	// get traction length and area
	double dl = sqrt(DotVectors2D(tang,tang));
	if(theCrack!=NULL)
		area = fmobj->IsAxisymmetric() ? cp.x*dl : theCrack->GetThickness()*dl ;
	
	// unit vector and length and get normal vector
	ScaleVector(tang,1./dl);
	*norm = MakeVector(-tang->y,tang->x,0.);
	
	// get nCod and tCod
	nCod = DotVectors(&cod,norm);
	tCod = DotVectors(&cod,tang);

	// return area (if needed)
	return area;
}

// update tractions
void CrackSegment::UpdateTractions(CrackHeader *theCrack)
{
	// exit if no traction law (matnum=0 or less return <0)
	if(MatID()<0) return;
	
	// get normal and tangential unit vectors and length
	Vector n,t;
	double nCod,tCod;
	double area = GetNormalAndTangent(theCrack,&n,&t,nCod,tCod);

	// will eventually call a traction law material and get total force
	// or force per radian if axisymmetric
	TractionLaw *theLaw=(TractionLaw *)theMaterials[MatID()];
	theLaw->CrackTractionLaw(this,nCod,tCod,&n,&t,area);
}

// get surface traction on this crack particle or zero (if no traction law)
// return true if set traction, or false is closed
bool CrackSegment::SurfaceTraction(CrackHeader *theCrack,Vector *ntract)
{
    // get normal and tangential unit vectors and length
    Vector n,t;
    double nCod,tCod;
    double area = GetNormalAndTangent(theCrack,&n,&t,nCod,tCod);
    
    // false if closed *may prefer small minimum here?)
    if(nCod<0.) return false;
    
    // zero if not using a traction law
    if(MatID()<0)
    {   ZeroVector(ntract);
        return true;
    }
    
    *ntract = tract;
    ScaleVector(ntract,-1./area);
    return true;
}


// Calculate energy in the traction law for this segment if in traction
// fullEnergy true gets total energy at location of this segment or nearest traction law tip
// fullEnergy false gets recoverable energy at the traction law tip
// Only used in J integral calculations for released and bridged energy, so return energy in N/mm
// If tipSegment is not NULL and this segment is not in the cohesive zone, it is set
//    to the crack segment at the beginning of the cohesive zone (or to NULL if no cohesive zone),
//    otherwise it is not changed
double CrackSegment::TractionEnergy(Vector *crossPt,int crkTipIdx,bool fullEnergy,CrackSegment **tipSegment)
{
	// if no traction law, scan toward crack tip looking for one
	if(MatID()<0)
	{	CrackSegment *closerSeg=this;
		while(TRUE)
		{	closerSeg = (crkTipIdx==START_OF_CRACK) ? closerSeg->prevSeg : closerSeg->nextSeg ;
			if(closerSeg==NULL) break;		// no traction law found near the crack tip
			
			// if find traction law, get crack tip energy (no need to interpolate)
			if(closerSeg->MatID()>=0)
			{	if(tipSegment!=NULL) *tipSegment = closerSeg;
				return closerSeg->SegmentTractionEnergy(fullEnergy);
			}
		}
		
		// no traction law on this crack; should return energy due to friction if crack surfaces in contact
		if(tipSegment!=NULL) *tipSegment = NULL;
		return 0.;
	}
	else if(!fullEnergy)
	{	// for released energy, scan to tip of this traction law
		CrackSegment *fartherSeg=this;
		while(TRUE)
		{	CrackSegment *closerSeg=fartherSeg;
			fartherSeg = (crkTipIdx==START_OF_CRACK) ? closerSeg->nextSeg : closerSeg->prevSeg ;
			if(fartherSeg==NULL) break;		// no traction law tip found (entire crack is traction law, may be non physical result)
			
			// if find end of traction law, get traction law tip recoverable energy from previous segment in the traction law
			if(fartherSeg->MatID()<0)
				return closerSeg->SegmentTractionEnergy(false);
		}
		
		// no traction law tip was found
		return 0.;
	}
	
	// get this segment's energy
	double energy=SegmentTractionEnergy(fullEnergy);
	
	// check other side of the cross point and interpolate if has energy too
	CrackSegment *fartherSeg = (crkTipIdx==START_OF_CRACK) ? nextSeg : prevSeg;
	if(fartherSeg!=NULL && crossPt!=NULL)
	{	if(fartherSeg->MatID()>=0)
		{	double energy2=fartherSeg->SegmentTractionEnergy(fullEnergy);
	
			// get segment line
			double dx=fartherSeg->cp.x-cp.x;
			double dy=fartherSeg->cp.y-cp.y;
			double dl=sqrt(dx*dx+dy*dy);
			
			// fraction of distance to the cross point
			double dcpx=crossPt->x-cp.x;
			double dcpy=crossPt->y-cp.y;
			double fract=sqrt(dcpx*dcpx+dcpy*dcpy)/dl;
			
			// interpolate energy
			energy=energy+fract*(energy2-energy);
		}
	}
	
	// return interpolated result
	return energy;
}

// Calculate potentialEnergy (if true) or dissipated energy (if false) for segment
// If not at crack tip, scan to one at the crack tip
// Used in J integral calculations so return energy in N/mm
double CrackSegment::SegmentTractionEnergy(bool getPotentialEnergy)
{
	// call on traction law material for this segment
	TractionLaw *theLaw=(TractionLaw *)theMaterials[MatID()];

	if(getPotentialEnergy)
	{	// get normal and tangential unit vector and length
		Vector n,t;
		double nCod,tCod;
		GetNormalAndTangent(NULL,&n,&t,nCod,tCod);
		return theLaw->CrackWorkEnergy(this,nCod,tCod);
	}
	else
	{	double GI,GII;
		theLaw->CrackDissipatedEnergy(this,GI,GII);
		return GI+GII;
	}
}

// fill archive with this object values
void CrackSegment::FillArchive(char *app,int segNum,bool threeD)
{
    char *appInit=app;
    
    // must have these defaults
    *(int *)app=planeInElem;
    app+=sizeof(int);
	
	// tipMatNum save as int, but in space for double (for backward compatibility)
	*(int *)app=tipMatnum;
	app+=sizeof(double);
    
	// segNum==0 at start of new crack
    // start new crack by archiving -1
	//  or continue crack by archiving -2
    if(segNum)
        *(short *)app=-2;
    else
        *(short *)app=-1;
    app+=sizeof(short);
	
	// 1-based material ID for traction law
	// warning: this value was undefined in files before tractions were added
	*(short *)app=matnum;
    app+=sizeof(short);
    
    // new defaults assumes archiver->crackOrder[ARCH_Defaults]=='Y'
    
    // position
    *(double *)app=cp.x;
    app+=sizeof(double);
    
    *(double *)app=cp.y;
    app+=sizeof(double);
    
    if(threeD)
    {   *(double *)app=cp.z;
        app+=sizeof(double);
    }
    
    // original position
    *(double *)app=orig.x;
    app+=sizeof(double);
    
    *(double *)app=orig.y;
    app+=sizeof(double);
    
    if(threeD)
    {   *(double *)app=orig.z;
        app+=sizeof(double);
    }
    
    // above the crack
    *(int *)app=surfInElem[0];
    app+=sizeof(int);
    
    *(double *)app=surf[0].x;
    app+=sizeof(double);
    
    *(double *)app=surf[0].y;
    app+=sizeof(double);
    
    if(threeD)
    {   *(double *)app=surf[0].z;
        app+=sizeof(double);
    }
    
    // below the crack
    *(int *)app=surfInElem[1];
    app+=sizeof(int);
    
    *(double *)app=surf[1].x;
    app+=sizeof(double);
    
    *(double *)app=surf[1].y;
    app+=sizeof(double);
    
    if(threeD)
    {   *(double *)app=surf[1].z;
        app+=sizeof(double);
    }
    
    // J integral
    if(archiver->CrackArchive(ARCH_JIntegral))
    {	*(double *)app=Jint.x*UnitsController::Scaling(1.e-3);		// current crack tip J1
        app+=sizeof(double);
		if(fmobj->propagate[0])
			*(double *)app=propagationJ*UnitsController::Scaling(1.e-3);		// actual energy released last time the crack grew
		else
			*(double *)app=Jint.y*UnitsController::Scaling(1.e-3);		// current crack tip J2
        app+=sizeof(double);
    }
    
    // Stress Intensity Factors (/sqrt(1000) for  units Pa sqrt(m))
    if(archiver->CrackArchive(ARCH_StressIntensity))
    {	*(double *)app=sif.x*UnitsController::Scaling(0.0316227766016838e-6);
        app+=sizeof(double);
        *(double *)app=sif.y*UnitsController::Scaling(0.0316227766016838e-6);
        app+=sizeof(double);
    }
	
	// Cummulative energy dissipated per unit length by this cohesize zone
	if(archiver->CrackArchive(ARCH_CZMDISP))
    {   if(czmdG.z>0.)
        {   *(int *)app=1;
            app+=sizeof(int);
            // Energy releaase rate in N/mm
            *(double *)app=czmdG.x*UnitsController::Scaling(1.e-6);
            app+=sizeof(double);
            *(double *)app=czmdG.y*UnitsController::Scaling(1.e-6);
            app+=sizeof(double);
        }
        else
        {   *(int *)app=0;
            app+=sizeof(int);
            *(double *)app=0.;
            app+=sizeof(double);
            *(double *)app=0.;
            app+=sizeof(double);
        }
	}
    
    // Traction 1 to 5
    int tractMat = MatID();
    char tractionChar = archiver->GetCrackOrderByte(ARCH_Traction15);
    if(tractionChar=='Y')
    {   *(double *)app = tractMat>=0? theMaterials[tractMat]->GetHistory(1,historyData) : 0. ;
        app+=sizeof(double);
    }
    else if(tractionChar!='N')
    {   if(tractionChar&0x01)
        {   *(double *)app = tractMat>=0? theMaterials[tractMat]->GetHistory(1,historyData) : 0. ;
            app+=sizeof(double);
        }
        if(tractionChar&0x02)
        {   *(double *)app = tractMat>=0? theMaterials[tractMat]->GetHistory(2,historyData) : 0. ;
            app+=sizeof(double);
        }
        if(tractionChar&0x04)
        {   *(double *)app = tractMat>=0? theMaterials[tractMat]->GetHistory(3,historyData) : 0. ;
            app+=sizeof(double);
        }
        if(tractionChar&0x08)
        {   *(double *)app = tractMat>=0? theMaterials[tractMat]->GetHistory(4,historyData) : 0. ;
            app+=sizeof(double);
        }
        if(tractionChar&0x10)
        {   *(double *)app = tractMat>=0? theMaterials[tractMat]->GetHistory(5,historyData) : 0. ;
            app+=sizeof(double);
        }
    }

    // Traction 6 to 10
    tractionChar = archiver->GetCrackOrderByte(ARCH_Traction610);
    if(tractionChar=='Y')
    {   *(double *)app = tractMat>=0? theMaterials[tractMat]->GetHistory(6,historyData) : 0. ;
        app+=sizeof(double);
    }
    else if(tractionChar!='N')
    {   if(tractionChar&0x01)
        {   *(double *)app = tractMat>=0? theMaterials[tractMat]->GetHistory(6,historyData) : 0. ;
            app+=sizeof(double);
        }
        if(tractionChar&0x02)
        {   *(double *)app = tractMat>=0? theMaterials[tractMat]->GetHistory(7,historyData) : 0. ;
            app+=sizeof(double);
        }
        if(tractionChar&0x04)
        {   *(double *)app = tractMat>=0? theMaterials[tractMat]->GetHistory(8,historyData) : 0. ;
            app+=sizeof(double);
        }
        if(tractionChar&0x08)
        {   *(double *)app = tractMat>=0? theMaterials[tractMat]->GetHistory(9,historyData) : 0. ;
            app+=sizeof(double);
        }
        if(tractionChar&0x10)
        {   *(double *)app = tractMat>=0? theMaterials[tractMat]->GetHistory(10,historyData) : 0. ;
            app+=sizeof(double);
        }
    }

    // reverse bytes if needed
    if(fmobj->GetReverseBytes())
    {	app=appInit;
    
        // must have these defaults
        app+=Reverse(app,sizeof(int));
		Reverse(app,sizeof(int));			// tipMatnum stored as int, but in space for double
        app+=sizeof(double);
        app+=Reverse(app,sizeof(short))+2;
        app+=Reverse(app,sizeof(double));	// position
        app+=Reverse(app,sizeof(double));
        app+=Reverse(app,sizeof(double));	// orig position
        app+=Reverse(app,sizeof(double));
        app+=Reverse(app,sizeof(int));		// above crack
        app+=Reverse(app,sizeof(double));
        app+=Reverse(app,sizeof(double));
        app+=Reverse(app,sizeof(int));		// below crack
        app+=Reverse(app,sizeof(double));
        app+=Reverse(app,sizeof(double));
        
		if(archiver->CrackArchive(ARCH_JIntegral))
        {   app+=Reverse(app,sizeof(double));
            app+=Reverse(app,sizeof(double));
        }
		if(archiver->CrackArchive(ARCH_StressIntensity))
        {   app+=Reverse(app,sizeof(double));
            app+=Reverse(app,sizeof(double));
        }
		if(archiver->CrackArchive(ARCH_CZMDISP))
		{   app+=Reverse(app,sizeof(int));
            app+=Reverse(app,sizeof(double));
            app+=Reverse(app,sizeof(double));
        }
        // traction history data 1 to 5
        tractionChar = archiver->GetCrackOrderByte(ARCH_Traction15);
        if(tractionChar=='Y')
            app+=Reverse(app,sizeof(double));
        else if(tractionChar!='N')
        {   if(tractionChar&0x01) app+=Reverse(app,sizeof(double));
            if(tractionChar&0x02) app+=Reverse(app,sizeof(double));
            if(tractionChar&0x04) app+=Reverse(app,sizeof(double));
            if(tractionChar&0x08) app+=Reverse(app,sizeof(double));
            if(tractionChar&0x10) app+=Reverse(app,sizeof(double));
        }
        // traction history data 6 to 10
        tractionChar = archiver->GetCrackOrderByte(ARCH_Traction610);
        if(tractionChar=='Y')
            app+=Reverse(app,sizeof(double));
        else if(tractionChar!='N')
        {   if(tractionChar&0x01) app+=Reverse(app,sizeof(double));
            if(tractionChar&0x02) app+=Reverse(app,sizeof(double));
            if(tractionChar&0x04) app+=Reverse(app,sizeof(double));
            if(tractionChar&0x08) app+=Reverse(app,sizeof(double));
            if(tractionChar&0x10) app+=Reverse(app,sizeof(double));
        }
    }
}

// Tell new crack tip to heat itself after the recent propagation
void CrackSegment::StartCrackTipHeating(double growth,double thickness)
{
	MaterialBase *tipMat=theMaterials[tipMatnum-1];
	double fractH = 1.0;								// fraction to heat (should be material property)
	
	// Will get wave speed of material. It would be better to have current wave speed, but
	// we don't have material point pointer to known and the crack tip will have particles
	// in different states
	double adot = tipMat->WaveSpeed(false,NULL);		// in mm/sec
	
	// set up rate and times (making sure proper number of steps)
	int nsteps = (int)(growth/(adot*timestep));
	if(nsteps<1) nsteps=1;
	
	// get heat rate E divided up among nsteps t get Pwr = E/T
	heatRate = fractH*Jint.x*thickness*growth/(timestep*(double)nsteps);
	
	// stop heating at this time
	// discrete version of heatEndTime=mtime + growth/adot;
	heatEndTime = mtime + ((double)nsteps+.2)*timestep;
	
	heating = true;
}

// return heating rate if any in nW/g
double CrackSegment::HeatRate(void)
{
	// return zero if off or if done
	if(!heating) return 0.;
	if(mtime>heatEndTime)
	{	heating = false;
		heatRate = 0.;
	}
	
	// return the heat rate
	return heatRate;
}

// return traction force (times a shape function) in uN
Vector CrackSegment::FTract(double fni)
{	Vector fout;
	fout.x=fni*tract.x;
	fout.y=fni*tract.y;
	fout.z=0.;
	return fout;
}

// calculate tractions on one side of crack for this segment
// throws std::bad_alloc
void CrackSegment::FindCrackTipMaterial(int currentNum)
{
	// if only one active material, it cannot change
	if(numActiveMaterials<=1) return;
	
	// skip if constant material
	if(theMaterials[currentNum-1]->KeepsCrackTip()) return;
	
	int i,iel;
#ifdef CONST_ARRAYS
	double fn[MAX_SHAPE_NODES];
	int nds[MAX_SHAPE_NODES];
#else
	double fn[maxShapeNodes];
	int nds[maxShapeNodes];
#endif

	// array to collect weights
	double *matWeight = new double[numActiveMaterials];
	for(i=0;i<numActiveMaterials;i++) matWeight[i] = 0.;
	
	// get shape functions to extrapolate to the particle
	iel=planeElemID();
	theElements[iel]->GetShapeFunctionsForCracks(fn,nds,cp);
	int numnds = nds[0];
	
	// extrapolate to particle
	for(i=1;i<=numnds;i++)
	{	nd[nds[i]]->AddMatWeights(fn[i],matWeight);
	}
	
	// find largest one
	int tipMat = -1;
	double tipWeight=0.;
	for(i=0;i<numActiveMaterials;i++)
	{	if(matWeight[i]>tipWeight)
		{	tipMat = i;
			tipWeight = matWeight[i];
		}
	}
	
	// free memory
	delete [] matWeight;
	
	// if found one, use it
	if(tipMat>=0)
		tipMatnum = MaterialBase::GetActiveMatID(tipMat)+1;

}

#pragma mark PREVENT PLANE CROSSES

// After crack plane moves, verify it has not passed either surface.
// If it has, slide that surface along COD to crack plane
// See JAN-OSU-4, pg 76 for algorithm. final equation on pg 82
int CrackSegment::CheckSurfaces(void)
{
	// check top surface
	if(prevSeg==NULL)
	{	// first segment only
		if(CrackHeader::Triangle(surf[0].x,surf[0].y,cp.x,cp.y,nextSeg->cp.x,nextSeg->cp.y)<0.)
		{	MoveToPlane(ABOVE_CRACK,nextSeg->cp.x-cp.x,nextSeg->cp.y-cp.y,false,1.);
			if(!FindElement(ABOVE_CRACK)) return false;
		}
	}
	else if(nextSeg==NULL)
	{	// last segment only
		if(CrackHeader::Triangle(surf[0].x,surf[0].y,prevSeg->cp.x,prevSeg->cp.y,cp.x,cp.y)<0.)
		{	MoveToPlane(ABOVE_CRACK,prevSeg->cp.x-cp.x,prevSeg->cp.y-cp.y,false,-1.);
			if(!FindElement(ABOVE_CRACK)) return false;
		}
	}
	else
	{	// internal segments, check each until path intersects the crack segment
		bool moved=false;
		if(CrackHeader::Triangle(surf[0].x,surf[0].y,prevSeg->cp.x,prevSeg->cp.y,cp.x,cp.y)<0.)
		{	moved=MoveToPlane(ABOVE_CRACK,prevSeg->cp.x-cp.x,prevSeg->cp.y-cp.y,true,-1.);
			if(moved)
			{	if(!FindElement(ABOVE_CRACK)) return false;
			}
		}
		if(!moved)
		{	if(CrackHeader::Triangle(surf[0].x,surf[0].y,cp.x,cp.y,nextSeg->cp.x,nextSeg->cp.y)<0.)
			{	if(MoveToPlane(ABOVE_CRACK,nextSeg->cp.x-cp.x,nextSeg->cp.y-cp.y,false,1.))
				{	if(!FindElement(ABOVE_CRACK)) return false;
				}
			}
		}
	}
	
	// check bottom surface
	if(prevSeg==NULL)
	{	// first segment only
		if(CrackHeader::Triangle(surf[1].x,surf[1].y,cp.x,cp.y,nextSeg->cp.x,nextSeg->cp.y)>0.)
		{	MoveToPlane(BELOW_CRACK,nextSeg->cp.x-cp.x,nextSeg->cp.y-cp.y,false,-1.);
			if(!FindElement(BELOW_CRACK)) return false;
		}
	}
	else if(nextSeg==NULL)
	{	// last segment only
		if(CrackHeader::Triangle(surf[1].x,surf[1].y,prevSeg->cp.x,prevSeg->cp.y,cp.x,cp.y)>0.)
		{	MoveToPlane(BELOW_CRACK,prevSeg->cp.x-cp.x,prevSeg->cp.y-cp.y,false,1.);
			if(!FindElement(BELOW_CRACK)) return false;
		}
	}
	else
	{	// internal segments
		bool moved=false;
		if(CrackHeader::Triangle(surf[1].x,surf[1].y,prevSeg->cp.x,prevSeg->cp.y,cp.x,cp.y)>0.)
		{	moved=MoveToPlane(BELOW_CRACK,prevSeg->cp.x-cp.x,prevSeg->cp.y-cp.y,true,1.);
			if(moved)
			{	if(!FindElement(BELOW_CRACK)) return false;
			}
		}
		if(!moved)
		{	if(CrackHeader::Triangle(surf[1].x,surf[1].y,cp.x,cp.y,nextSeg->cp.x,nextSeg->cp.y)>0.)
			{	if(MoveToPlane(BELOW_CRACK,nextSeg->cp.x-cp.x,nextSeg->cp.y-cp.y,false,-1.))
				{	if(!FindElement(BELOW_CRACK)) return false;
				}
			}
		}
	}
	
	return true;
}

// Move surface side (ABOVE_CRACK or BELOW_CRACK) to plane
// See JAN-OSU-4, pg 133-134
// thereIsAnotherSegement only true on first of two checks for internal segments
// vector dir*(-dyp,dxp) should point from surface to back across the crack normal to segment.
// side = ABOVE_CRACK (1) or BELOW)_CRACK (2)
bool CrackSegment::MoveToPlane(int side,double dxp,double dyp,bool thereIsAnotherSegement,double dir)
{
	int j=side-1;								// index to surface position
	double dxs=surf[j].x-cp.x,dys=surf[j].y-cp.y;		// vector crack particle to surface particle = xs-x
	double segLength=sqrt(dxp*dxp+dyp*dyp);		// segment length
	
	// distance to crack particle (relative to segment length)
	// if far away, then ignore it. Far away means might have misintrepreted the side
	double cod=sqrt(dxs*dxs+dys*dys)/segLength;
	if(cod>1.) return !thereIsAnotherSegement;

	// distance crack particle to intersection plane (relative to segment length)
	double t=(dxs*dxp+dys*dyp)/segLength;
	
	// if less than 0 and at first of two internal segments, return to try the next segment instead
	if(t<0. && thereIsAnotherSegement) return false;
	
	// distance surface to crack plane (relative to segment length)
	double n=(dxs*dyp-dys*dxp)*dir/segLength;
	
	// if pretty close or negative, then do not move to plane and hope velocity fields will resolve on their own
	// or if first of two internal segments, try the other one
	if(n<1.e-5) return !thereIsAnotherSegement;
	
	// restrict terminal segments
	if(prevSeg==NULL || nextSeg==NULL)
    {   if(t<0.) t=0.;
	}
	
	// restrict to intersect within this segment
	if(t>1.)
		t=1.;
	else if(t<-1.)
		t=-1.;
	
	// move to crack plane and a little more in normal direction
	surf[j].x=cp.x+t*dxp-1.0e-12*dir*dyp;
	surf[j].y=cp.y+t*dyp+1.0e-12*dir*dxp;
	
	return true;
}

#pragma mark HIERACHICAL CRACKS

// Create extents when segment is created or changed
// This is the first segment if isFirstSeg is TRUE and it is the last
//  if nextSeg->nextSeg is NULL
void CrackSegment::CreateSegmentExtents(bool isFirstSeg)
{
	if(nextSeg==NULL) return;           // end of the crack so no extents needed
    
    // fetch endpoints of line from this crack particle to the next
    double x1 = cp.x;
    double y1 = cp.y;
    double x2 = nextSeg->cp.x;
    double y2 = nextSeg->cp.y;
    
    // check for exterior crack on either end of the crack
    if(isFirstSeg && tipMatnum==EXTERIOR_CRACK)
    {   x1 = 5.*x1-4.*x2;		// x1-4*(x2-x1)
        y1 = 5.*y1-4.*y2;		// y1-4*(y2-y1)
    }
    if(nextSeg->nextSeg==NULL)
    {   if(nextSeg->tipMatnum==EXTERIOR_CRACK)
        {   x2=5.*x2-4.*x1;		// xn+4*(xn-x(nm1))
            y2=5.*y2-4.*y1;		// yn+4*(yn-y(nm1))
        }
    }
    
    // Find extents for bounding octagon
    
    // Pi = (1,0)
    if(x1>x2)
    {   cnear[0] = x2;        // min
        cfar[0] = x1;         // max
    }
    else
    {   cnear[0] = x1;        // min
        cfar[0] = x2;         // max
    }
    
    // Pi = (0,1)
    if(y1>y2)
    {   cnear[1] = y2;        // min
        cfar[1] = y1;         // max
    }
    else
    {   cnear[1] = y1;        // min
        cfar[1] = y2;         // max
    }
    
    // Pi = (1,1)
    double sum1=x1+y1,sum2=x2+y2;
    if(sum1>sum2)
    {   cnear[2] = sum2;        // min
        cfar[2] = sum1;         // max
    }
    else
    {   cnear[2] = sum1;        // min
        cfar[2] = sum2;         // max
    }
    
    // Pi = (1,-1)
    double diff1=x1-y1,diff2=x2-y2;
    if(diff1>diff2)
    {   cnear[3] = diff2;        // min
        cfar[3] = diff1;         // max
    }
    else
    {   cnear[3] = diff1;        // min
        cfar[3] = diff2;         // max
    }
}

#pragma mark ACCESSORS

// move above or below position slightly in the direction of the normal
// if it has already moved, then do not bother
// side = ABOVE_CRACK (1) or BELOW)_CRACK (2)
Vector CrackSegment::SlightlyMovedIfNotMovedYet(int side)
{
	Vector moved;
	moved = surf[side-1];
	if(cp.x!=moved.x && cp.y!=moved.y) return moved;
	
    // get vector normal to crack from above to below
	// v = (dx,dy) is segment (to next or from previous), then n = (dy,-dx)
	double dx,dy;
	if(nextSeg!=NULL)
	{	dy=cp.x-nextSeg->cp.x;		// -dx
		dx=nextSeg->cp.y-cp.y;		// dy
	}
	else
	{	dy=prevSeg->cp.x-cp.x;		// -dx
		dx=cp.y-prevSeg->cp.y;		// dy
	}
	
    // move in direction of this side
	if(side==ABOVE_CRACK)
	{	moved.x-=dx*1e-8;
		moved.y-=dy*1e-8;
	}
	else
	{	moved.x+=dx*1e-8;
		moved.y+=dy*1e-8;
	}
	
	return moved;
}
// material ID (convert to zero based)
int CrackSegment::MatID(void) { return matnum-1; }			// convert 1-based matnum to zero based for materials array
void CrackSegment::SetMatID(int newMat) { matnum=newMat; }			// input 1-based

// pointer to history data for use by traction laws
void CrackSegment::SetHistoryData(char *p)
{	if(historyData!=NULL) delete [] historyData;
	historyData=p;
}
char *CrackSegment::GetHistoryData(void) { return historyData; }

// zero based element ID to addess into elements array
// side = ABOVE_CRACK (1) or BELOW)_CRACK (2)
int CrackSegment::planeElemID(void) const { return planeInElem-1; }
int CrackSegment::surfaceElemID(int side) const { return surfInElem[side-1]-1; }

// store crack header (2D segment does not)
void CrackSegment::SetCrackHeader(CrackHeader *theHeader) { header = theHeader; }

