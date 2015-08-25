/********************************************************************************
    CrackSegment.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Apr 3 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Cracks/CrackSegment.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Elements/ElementBase.hpp"
#include "Materials/TractionLaw.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "System/ArchiveData.hpp"
#include "Cracks/CrackHeader.hpp"
#include "System/UnitsController.hpp"

extern char *app;

#pragma mark CrackSegment: Constructors and Destructors

// Constructors
CrackSegment::CrackSegment()
{
}

CrackSegment::CrackSegment(double xend,double yend,int tip,int matid)
{
    x=origx=surfx[0]=surfx[1]=xend;
    y=origy=surfy[0]=surfy[1]=yend;
    planeInElem=surfInElem[0]=surfInElem[1]=0;
    nextSeg=NULL;
	prevSeg=NULL;
    parent=NULL;
    
	ZeroVector(&Jint);
	ZeroVector(&sif);
	ZeroVector(&tract);
    tipMatnum=tip;
	propagationJ=0.;
    steadyState=STATIONARY;
	heating=FALSE;
	SetMatID(matid);
	historyData=NULL;
}

#pragma mark CrackSegment: Methods

// find current element (1 based) or return 0 if no element
int CrackSegment::FindElement(void)
{
    int i;
	Vector cpt;
	cpt.x=x;
	cpt.y=y;
	cpt.z=0.;
    
    // check current element
    if(planeInElem>0)
    {	if(theElements[planeInElem-1]->PtInElement(cpt))
            return planeInElem;
    }
    
    // check others
    planeInElem=0;
    for(i=0;i<nelems;i++)
    {	if(theElements[i]->PtInElement(cpt))
        {   planeInElem=i+1;
			if(theElements[i]->OnTheEdge()) planeInElem=0;
            break;
        }
    }
    return planeInElem;
}

// find current element (1 based) or return 0 if no element for crack surface
int CrackSegment::FindElement(short side)
{
    int i,j=side-1;
	Vector cpt;
	cpt.x=surfx[j];
	cpt.y=surfy[j];
	cpt.z=0.;
	
    // check current element
    if(surfInElem[j]>0)
    {	if(theElements[surfInElem[j]-1]->PtInElement(cpt))
            return surfInElem[j];
    }
	
    // check others
    surfInElem[j]=0;
    for(i=0;i<nelems;i++)
    {	if(theElements[i]->PtInElement(cpt))
        {   surfInElem[j]=i+1;
			if(theElements[i]->OnTheEdge()) surfInElem[j]=0;
            break;
        }
    }
    return surfInElem[j];
}

// Reset crack plane position from surfaces (2D) (in mm)
// only used when contact.GetMoveOnlySurfaces() is true and thus get crack position
//		from previous movement of the surfaces
void CrackSegment::MovePosition(void)
{	
	x += dxPlane;
	y += dyPlane;
    //x=(surfx[0]+surfx[1])/2.;
    //y=(surfy[0]+surfy[1])/2.;
}

// Reset crack plane position from surfaces (2D) (in mm)
void CrackSegment::MovePositionToMidpoint(void)
{
    x=(surfx[0]+surfx[1])/2.;
    y=(surfy[0]+surfy[1])/2.;
}

// Move crack plane position (2D) (in mm)
// only used when contact.GetMoveOnlySurfaces() is FALSE and thus need to move crack plane
void CrackSegment::MovePosition(double xpt,double ypt)
{	
    x+=xpt;
    y+=ypt;
}

// Move a surface position (2D) (in mm) - must move ABOVE_CRACK and then BELOW_CRACK
bool CrackSegment::MoveSurfacePosition(short side,double xpt,double ypt,bool hasNodes)
{
    short j=side-1;
	bool movedOther=FALSE;
	
	if(side==ABOVE_CRACK)
	{	// above crack is first
		if(hasNodes)
		{	surfx[j] += xpt;			// move
			surfy[j] += ypt;
			dxPlane = xpt;				// save until below is done next
			dyPlane = ypt;
		}
		else
		{	dxPlane = 0.;
			dyPlane = 0.;
		}
		hadAboveNodes = hasNodes;
	}
	else
	{	// below crack is second
		// ... if has nodes than get average movement or just this one if top did not move
		//			also move top along if it had no nodes
		// ... else if no nodes, move along with the top field (if it had one)
		if(hasNodes)
		{	surfx[j] += xpt;			// move
			surfy[j] += ypt;
			
			if(hadAboveNodes)
			{	dxPlane = 0.5*(dxPlane+xpt);
				dyPlane = 0.5*(dyPlane+ypt);
			}
			else
			{	dxPlane = xpt;					// only had nodes below the crack
				dyPlane = ypt;
				surfx[ABOVE_CRACK-1] += xpt;	// move above by (xpt,ypt) too
				surfy[ABOVE_CRACK-1] += ypt;
				movedOther = true;				// in case above moved elements now (this tells to check on return)
			}
		}
		else if(hadAboveNodes)
		{	surfx[j] += dxPlane;				// only had nodes above the crack, move both and th plane by (dxPlane,dyPlane)
			surfy[j] += dyPlane;
		}
	}
	
	return movedOther;
}

// calculate tractions due to this segment
void CrackSegment::AddTractionForceSeg(CrackHeader *theCrack)
{	// exit if no traction law
	if(MatID()<0) return;
	
	// The first call will add T to each node and find fnorm, which is
	//    sum of f_i shape functions for all nodes, but to account for partia;
	//    nodes, it should have added (f_i/fnorm) to each node.
	// Thus if fnorm is not close to 1 (or 0, which means no nodes gound), call
	//    again scaling by 1/fnorm - 1. The first term gets the correct f_i/fnorm
	//    and the -1 subtracts out the first call
	// The net result will be application of T regardless of number of nodes
	double fnorm=AddTractionForceSegSide(theCrack,ABOVE_CRACK,(double)1.);
	if(fnorm>0. && fnorm<0.999)
		AddTractionForceSegSide(theCrack,ABOVE_CRACK,(1./fnorm-1.));
	
	// Repeat above calculations for the below the crack field, with -T
	fnorm=AddTractionForceSegSide(theCrack,BELOW_CRACK,(double)-1.);
	if(fnorm>0. && fnorm<0.999)
		AddTractionForceSegSide(theCrack,ABOVE_CRACK,(1.-1./fnorm));
}

// calculate tractions on one side of crack for this segment
// add forces to material velocity fields on one side of the crack
double CrackSegment::AddTractionForceSegSide(CrackHeader *theCrack,int side,double sign)
{
	int numnds,nds[maxShapeNodes];
    double fn[maxShapeNodes];
    short vfld;
	double fnorm = 0.;
	NodalPoint *ndi;
	int cnum=theCrack->GetNumber();

	Vector cspos = MakeVector(surfx[side-1], surfy[side-1], 0.);
	const ElementBase *elref = theElements[surfInElem[side-1]];
	elref->GetShapeFunctionsForCracks(&numnds,fn,nds,&cspos);
	
	// loop over all nodes seen by this crack surface particle
	for(int i=1;i<=numnds;i++)
	{	// Get velocity field to use
		ndi = nd[nds[i]];
		vfld = ndi->GetFieldForSurfaceParticle(side,cnum,this);
		
		// if has particles (and they see cracks), add force and track normalization
		if(vfld>=0)
		{	ndi->AddFtotSpreadTask3(vfld,FTract(sign*fn[i]));
			fnorm += fn[i];
		}
	}
	
	// return amount used
	return fnorm;
}

// get normalized vector tangent to crack path and the length of the path associated with this particle
// vector points in direction of the crack. The corresponding normal (-t.y,t.x) is from below to above
// (better to use splines)
Vector CrackSegment::GetTangential(double *length)
{
	Vector tang;
	double dl;
	
	if(prevSeg==NULL)
	{	tang.x=nextSeg->x-x;
		tang.y=nextSeg->y-y;
	}
	else if(nextSeg==NULL)
	{	tang.x=x-prevSeg->x;
		tang.y=y-prevSeg->y;
	}
	else
	{	tang.x=(nextSeg->x-prevSeg->x)/2.;
		tang.y=(nextSeg->y-prevSeg->y)/2.;
	}
	dl=sqrt(tang.x*tang.x+tang.y*tang.y);
	
	// unit vector and length
	tang.x/=dl;
	tang.y/=dl;
	tang.z=0.;
	*length=dl;
	
	return tang;
}

// update tractions
void CrackSegment::UpdateTractions(CrackHeader *theCrack)
{	
	// exit if no traction law (matnum=0 or less return <0)
	if(MatID()<0) return;
	
	// get tangential unit vector and length
	double dl;
	Vector t = GetTangential(&dl);
	
	// cutting assume vector along x axis, use length from above
	if(fmobj->dflag[0]==4)
	{	t.x=-1.;
		t.y=0.;
		t.z=0.;
	}
	
	// get normal and tangential COD components
	double codx=surfx[0]-surfx[1];			// above crack - below crack
	double cody=surfy[0]-surfy[1];
	double nCod=-t.y*codx+t.x*cody;			// absolute normal cod (delta.n where n from below to above)
	double tCod=t.x*codx+t.y*cody;			// absolute tangential cod (delta.t)
	
	// will eventually call a traction law material and get total force
	// or force per radian if axisymmetric
	double area = fmobj->IsAxisymmetric() ? x*dl : theCrack->GetThickness()*dl ;
	TractionLaw *theLaw=(TractionLaw *)theMaterials[MatID()];
	theLaw->CrackTractionLaw(this,nCod,tCod,t.x,t.y,area);

}

// Calculate energy in the traction law for this segment if in traction
// fullEnergy true gets total energy at location of this segment or nearest traction law tip
// full Energy false gets recoverable energy at the traction law tip
// Only used in J integral calculations for released and bridged energy, so return energy in N/mm
// If tipSegment is not NULL and this segment is not in the cohesive zone, it is set
//    to the crack segment and the beginning of the cohesive zone (or to NULL if no cohesive zone),
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
	if(fartherSeg!=NULL)
	{	if(fartherSeg->MatID()>=0)
		{	double energy2=fartherSeg->SegmentTractionEnergy(fullEnergy);
	
			// get segment line
			double dx=fartherSeg->x-x;
			double dy=fartherSeg->y-y;
			double dl=sqrt(dx*dx+dy*dy);
			
			// fraction of distance to the cross point
			double dcpx=crossPt->x-x;
			double dcpy=crossPt->y-y;
			double fract=sqrt(dcpx*dcpx+dcpy*dcpy)/dl;
			
			// interpolate energy
			energy=energy+fract*(energy2-energy);
		}
	}
	
	// return interpolated result
	return energy;
}

// Calculate energy in the traction law for this segment, assumed to be a traction material
// Used in J integral calculations so return energy in N/mm
double CrackSegment::SegmentTractionEnergy(bool fullEnergy)
{	
	// get tangential unit vector and length
	double dl;
	Vector t=GetTangential(&dl);
	
	// get normal and tangential COD components
	double codx=surfx[0]-surfx[1];			// above crack - below crack
	double cody=surfy[0]-surfy[1];
	double nCod=-t.y*codx+t.x*cody;			// absolute normal cod
	double tCod=t.x*codx+t.y*cody;			// absolute tangential cod
	
	// call on traction law material for this segment
	TractionLaw *theLaw=(TractionLaw *)theMaterials[MatID()];
	return theLaw->CrackTractionEnergy(this,nCod,tCod,fullEnergy);
}

// fill archive with this object values
void CrackSegment::FillArchive(char *app,int segNum)
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
    *(double *)app=x;
    app+=sizeof(double);
    
    *(double *)app=y;
    app+=sizeof(double);
    
    // original position
    *(double *)app=origx;
    app+=sizeof(double);
    
    *(double *)app=origy;
    app+=sizeof(double);
    
    // above the crack
    *(int *)app=surfInElem[0];
    app+=sizeof(int);
    
    *(double *)app=surfx[0];
    app+=sizeof(double);
    
    *(double *)app=surfy[0];
    app+=sizeof(double);
    
    // below the crack
    *(int *)app=surfInElem[1];
    app+=sizeof(int);
    
    *(double *)app=surfx[1];
    app+=sizeof(double);
    
    *(double *)app=surfy[1];
    app+=sizeof(double);
    
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
	
	// Energy Balance Results (no longer used)
	if(archiver->CrackArchive(ARCH_BalanceResults))
	{   *(int *)app=0;
		app+=sizeof(double);
        *(double *)app=0.;
        app+=sizeof(double);
        *(double *)app=0.;
        app+=sizeof(double);
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
		if(archiver->CrackArchive(ARCH_BalanceResults))
		{   app+=Reverse(app,sizeof(int));
            app+=Reverse(app,sizeof(double));
            app+=Reverse(app,sizeof(double));
        }
    }
}

// Tell new crack tip to heat itself after the recent propagation
void CrackSegment::StartCrackTipHeating(double growth,double thickness)
{
	MaterialBase *tipMat=theMaterials[tipMatnum-1];
	double fractH = 1.0;								// fraction to heat (should be material property)
	double adot = tipMat->WaveSpeed(FALSE,NULL);		// in mm/sec
	
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
void CrackSegment::FindCrackTipMaterial(int currentNum)
{
	// if only one active material, it cannot change
	if(numActiveMaterials<=1) return;
	
	// skip if constant material
	if(theMaterials[currentNum-1]->KeepsCrackTip()) return;
	
	Vector cspos;
	int i,iel,numnds,nds[maxShapeNodes];
    double fn[maxShapeNodes];
	
	// array to collect weights
	double *matWeight=(double *)malloc(sizeof(double)*numActiveMaterials);
	for(i=0;i<numActiveMaterials;i++) matWeight[i] = 0.;
	
	// get shape functions
	cspos.x=x;
	cspos.y=y;
	iel=planeInElem;
	theElements[iel]->GetShapeFunctionsForCracks(&numnds,fn,nds,&cspos);
	
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
		if(CrackHeader::Triangle(surfx[0],surfy[0],x,y,nextSeg->x,nextSeg->y)<0.)
		{	MoveToPlane(ABOVE_CRACK,nextSeg->x-x,nextSeg->y-y,false,1.);
			if(!FindElement(ABOVE_CRACK)) return false;
		}
	}
	else if(nextSeg==NULL)
	{	// last segment only
		if(CrackHeader::Triangle(surfx[0],surfy[0],prevSeg->x,prevSeg->y,x,y)<0.)
		{	MoveToPlane(ABOVE_CRACK,prevSeg->x-x,prevSeg->y-y,false,-1.);
			if(!FindElement(ABOVE_CRACK)) return false;
		}
	}
	else
	{	// internal segments, check each until path intersects the crack segment
		bool moved=false;
		if(CrackHeader::Triangle(surfx[0],surfy[0],prevSeg->x,prevSeg->y,x,y)<0.)
		{	moved=MoveToPlane(ABOVE_CRACK,prevSeg->x-x,prevSeg->y-y,true,-1.);
			if(moved)
			{	if(!FindElement(ABOVE_CRACK)) return false;
			}
		}
		if(!moved)
		{	if(CrackHeader::Triangle(surfx[0],surfy[0],x,y,nextSeg->x,nextSeg->y)<0.)
			{	if(MoveToPlane(ABOVE_CRACK,nextSeg->x-x,nextSeg->y-y,false,1.))
				{	if(!FindElement(ABOVE_CRACK)) return false;
				}
			}
		}
	}
	
	// check bottom surface
	if(prevSeg==NULL)
	{	// first segment only
		if(CrackHeader::Triangle(surfx[1],surfy[1],x,y,nextSeg->x,nextSeg->y)>0.)
		{	MoveToPlane(BELOW_CRACK,nextSeg->x-x,nextSeg->y-y,false,-1.);
			if(!FindElement(BELOW_CRACK)) return false;
		}
	}
	else if(nextSeg==NULL)
	{	// last segment only
		if(CrackHeader::Triangle(surfx[1],surfy[1],prevSeg->x,prevSeg->y,x,y)>0.)
		{	MoveToPlane(BELOW_CRACK,prevSeg->x-x,prevSeg->y-y,false,1.);
			if(!FindElement(BELOW_CRACK)) return false;
		}
	}
	else
	{	// internal segments
		bool moved=false;
		if(CrackHeader::Triangle(surfx[1],surfy[1],prevSeg->x,prevSeg->y,x,y)>0.)
		{	moved=MoveToPlane(BELOW_CRACK,prevSeg->x-x,prevSeg->y-y,true,1.);
			if(moved)
			{	if(!FindElement(BELOW_CRACK)) return false;
			}
		}
		if(!moved)
		{	if(CrackHeader::Triangle(surfx[1],surfy[1],x,y,nextSeg->x,nextSeg->y)>0.)
			{	if(MoveToPlane(BELOW_CRACK,nextSeg->x-x,nextSeg->y-y,false,-1.))
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
bool CrackSegment::MoveToPlane(int side,double dxp,double dyp,bool thereIsAnotherSegement,double dir)
{	
	int j=side-1;								// index to surface position
	double dxs=surfx[j]-x,dys=surfy[j]-y;		// vector crack particle to surface particle = xs-x
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
	surfx[j]=x+t*dxp-1.0e-12*dir*dyp;
	surfy[j]=y+t*dyp+1.0e-12*dir*dxp;
	
	return true;
}

// For cutting, collapse the crack if opened enough
// cods are hard coded and should be greater than traction law max cod
// Also won't collapse if this segment or a neighbor has traction law
// return tru unless collapsed position is out of the grid
bool CrackSegment::CollapseSurfaces(void)
{
	// make sure no traction on this segment or its neighbors
	if(MatID()>=0) return true;
	if(nextSeg!=NULL)
	{	if(nextSeg->MatID()>=0) return true;
	}
	if(prevSeg!=NULL)
	{	if(prevSeg->MatID()>=0) return true;
	}
	
	// check x and y opening
	double codx=fabs(surfx[0]-surfx[1]);
	double cody=fabs(surfy[0]-surfy[1]);
	if(codx>1.02 || cody>1.02)
	{	// OK to collapse
		x=(surfx[0]+surfx[1])/2.;
		y=(surfy[0]+surfy[1])/2.;
		if(!FindElement()) return false;
		surfx[0]=surfx[1]=x;
		surfy[0]=surfy[1]=y;
		surfInElem[0]=surfInElem[1]=planeInElem;
		SetMatID(-1);								// marks as collapsed
	}
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
    double x1 = x;
    double y1 = y;
    double x2 = nextSeg->x;
    double y2 = nextSeg->y;
    
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
Vector CrackSegment::SlightlyMovedIfNotMovedYet(int side)
{
	Vector moved;
	moved.x=surfx[side-1];
	moved.y=surfy[side-1];
	moved.z=0;
	if(x!=moved.x && y!=moved.y) return moved;
	
    // get vector normal to crack from above to below v = (dx,dy)
	double dx,dy;
	if(nextSeg!=NULL)
	{	dy=x-nextSeg->x;		// -¶x
		dx=nextSeg->y-y;		// ¶y
	}
	else
	{	dy=prevSeg->x-x;		// -¶x
		dx=y-prevSeg->y;		// ¶y
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


// area of triangle with given point and vector in forward crack direction along one surface
// assumes either prevSeg or nextSeg != NULL
double CrackSegment::ForwardArea(double xpt,double ypt,int side)
{	int js=side-1;
	return (prevSeg!=NULL) ?
		CrackHeader::Triangle(xpt,ypt,prevSeg->surfx[js],prevSeg->surfy[js],surfx[js],surfy[js]) :
		CrackHeader::Triangle(xpt,ypt,surfx[js],surfy[js],nextSeg->surfx[js],nextSeg->surfy[js]) ;
}

// material ID (convert to zero based)
int CrackSegment::MatID(void) { return matnum-1; }			// convert 1-based matnum to zero based for materials array
void CrackSegment::SetMatID(int newMat) { matnum=newMat; }	// input 1-based

// pointer to history data for use by traction laws
void CrackSegment::SetHistoryData(char *p)
{	if(historyData!=NULL) free(historyData);
	historyData=p;
}
char *CrackSegment::GetHistoryData(void) { return historyData; }

