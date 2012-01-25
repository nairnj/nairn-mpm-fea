/********************************************************************************
    CrackSegment.cpp
    NairnMPM
    
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
	ZeroVector(&Jint);
	ZeroVector(&sif);
	ZeroVector(&tract);
    tipMatnum=tip;
	crackIncrements=0;
	release=absorb=propagationJ=0.;
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
// only used when contact.GetMoveOnlySurfaces() is TRUE and thus get crack position
//		from previous movement of the surfaces
// return true or false in bothSurfaces if both surfaces actually moved
void CrackSegment::MovePosition(void)
{	
	x+=dxPlane;
	y+=dyPlane;
    //x=(surfx[0]+surfx[1])/2.;
    //y=(surfy[0]+surfy[1])/2.;
}

// Move crack plane position (2D) (in mm)
// only used when contact.GetMoveOnlySurfaces() is FALSE and thus need to move crack plane
void CrackSegment::MovePosition(double xpt,double ypt)
{	
    x+=xpt;
    y+=ypt;
}

// Move a surface position (2D) (in mm) - must move ABOVE_CRACK and then BELOW_CRACK
bool CrackSegment::MoveSurfacePosition(short side,double xpt,double ypt,bool hasNodes,double surfaceMass)
{
    short j=side-1;
	bool movedOther=FALSE;
	
	if(side==ABOVE_CRACK)
	{	// above crack is first
		if(hasNodes)
		{	surfx[j]+=xpt;			// move
			surfy[j]+=ypt;
			dxPlane=xpt;			// save until below is done next
			dyPlane=ypt;
			aboveMass=surfaceMass;
		}
		else
		{	dxPlane=0.;
			dyPlane=0.;
			aboveMass=0.;
		}
		hadAboveNodes=hasNodes;
	}
	else
	{	// below crack is second
		// ... if has nodes than get average movement or just this one if top did not move
		//			also move top along if it had no nodes
		// ... else if no nodes, move along with the top field (if it had one)
		if(hasNodes)
		{	surfx[j]+=xpt;			// move
			surfy[j]+=ypt;
			
			if(hadAboveNodes)
			{	double sumMass=aboveMass+surfaceMass;
				dxPlane=(aboveMass*dxPlane+surfaceMass*xpt)/sumMass;	// had nodes above and below the crack, find the average
				dyPlane=(aboveMass*dyPlane+surfaceMass*ypt)/sumMass;
			}
			else
			{	dxPlane=xpt;					// internal point only had nodes below the crack
				dyPlane=ypt;
				surfx[ABOVE_CRACK-1]+=xpt;		// move above by (xpt,ypt) too
				surfy[ABOVE_CRACK-1]+=ypt;
				movedOther=TRUE;				// in case it moved elements too
			}
		}
		else if(hadAboveNodes)
		{	surfx[j]+=dxPlane;				// only had nodes above the crack, move both by (dxPlane,dyPlane)
			surfy[j]+=dyPlane;
		}
	}
	
	return movedOther;
}

// calculate tractions due to this segment
void CrackSegment::AddTractionFext(CrackHeader *theCrack)
{	// exit if no traction law
	if(MatID()<0) return;
	
	// the second call makes sure the full traction force is applied to both sides
	// of the crack. It may not be needed (or even be good?)
	double fnorm=AddTractionFextSide(theCrack,ABOVE_CRACK,(double)1.);
	if(fnorm>0. && fnorm<0.999)
		AddTractionFextSide(theCrack,ABOVE_CRACK,(1./fnorm-1.));
	
	fnorm=AddTractionFextSide(theCrack,BELOW_CRACK,(double)-1.);
	if(fnorm>0. && fnorm<0.999)
		AddTractionFextSide(theCrack,ABOVE_CRACK,(1.-1./fnorm));
}

// calculate tractions on one side of crack for this segment
double CrackSegment::AddTractionFextSide(CrackHeader *theCrack,int side,double sign)
{
	side--;			// convert to 0 or 1 for above and below
	Vector cspos,ndpos,norm;
	int i,iel,numnds,nds[MaxShapeNds];
    short vfld;
    double fn[MaxShapeNds];
	NodalPoint *ndi;
	double fnorm=0.;
	
	cspos.x=surfx[side];
	cspos.y=surfy[side];
	iel=surfInElem[side];
	theElements[iel]->GetShapeFunctions(&numnds,fn,nds,&cspos,&ndpos);
	for(i=1;i<=numnds;i++)
	{	ndi=nd[nds[i]];
		vfld=theCrack->CrackCross(cspos.x,cspos.y,ndi->x,ndi->y,&norm);
		if(vfld>NO_CRACK)
		{	// a crossing field - to use it, must find correct field and crack number in a velocity field
			// traction laws may not handle multiple cracks or interacting fields correctly, but try to do something
			//  1. Possible: [0], [1], [3], [0]&[3], [1]&[2], [0]&[1], [1]&[3], [0]&[1]&[2],
			//			[0]&[1]&[3], [1]&[2]&[3], and [0]&[1]&[2]&[3]
			//  2. Never occurs [2], [0]&[2], [2]&[3], [0]&[2]&[3]
			int cnum=theCrack->GetNumber();
			if(CrackVelocityField::ActiveNonrigidField(ndi->cvf[1]))
			{	// Node with: [1], [1]&[2], [0]&[1], [1]&[3], [0]&[1]&[2], [0]&[1]&[3], [1]&[2]&[3], or [0]&[1]&[2]&[3]
				if(vfld==ndi->cvf[1]->location(FIRST_CRACK) && cnum==ndi->cvf[1]->crackNumber(FIRST_CRACK))
				{	// this means field one is correct - single crack calcs will always get here if the field is active
					vfld=1;
				}
				else if(CrackVelocityField::ActiveNonrigidField(ndi->cvf[2]))
				{	// Node with: [1]&[2], [0]&[1]&[2], [1]&[2]&[3], or [0]&[1]&[2]&[3]
					if(vfld==ndi->cvf[2]->location(FIRST_CRACK) && cnum==ndi->cvf[2]->crackNumber(FIRST_CRACK))
					{	// this means there are two cracks here and this crack should add to field [2]
						vfld=2;
					}
					else if(CrackVelocityField::ActiveNonrigidField(ndi->cvf[3]))
					{	if((vfld==ndi->cvf[3]->location(FIRST_CRACK) && cnum==ndi->cvf[3]->crackNumber(FIRST_CRACK)) || 
							(vfld==ndi->cvf[3]->location(SECOND_CRACK) && cnum==ndi->cvf[3]->crackNumber(SECOND_CRACK)))
						{	// this means there are two cracks here and this crack should add to field [3]
							vfld=3;
						}
					}
				}
				else if(CrackVelocityField::ActiveNonrigidField(ndi->cvf[3]))
				{	// Node with: [1], [0]&[1], [1]&[3], or [0]&[1]&[3]
					if((vfld==ndi->cvf[3]->location(FIRST_CRACK) && cnum==ndi->cvf[3]->crackNumber(FIRST_CRACK)) || 
					   (vfld==ndi->cvf[3]->location(SECOND_CRACK) && cnum==ndi->cvf[3]->crackNumber(SECOND_CRACK)))
					{	// this means there are two cracks here and this crack should add to field [3]
						vfld=3;
					}
				}
				else
				{	// none are correct, sometimes [0] is correct in this case
					// Use [0], but might want a warning to be issued
					vfld=0;
				}
			}
			else if(CrackVelocityField::ActiveNonrigidField(ndi->cvf[3]))
			{	// Node with: [3], [0]&[3]
				if((vfld==ndi->cvf[3]->location(FIRST_CRACK) && cnum==ndi->cvf[3]->crackNumber(FIRST_CRACK)) || 
				   (vfld==ndi->cvf[3]->location(SECOND_CRACK) && cnum==ndi->cvf[3]->crackNumber(SECOND_CRACK)))
				{	// this means there are two cracks here and this crack should add to field [3]
					vfld=3;
				}
				else
				{	// none are correct, sometimes [0] is correct in this case
					// Use [0], but might want a warning to be issued
					vfld=0;
				}
			}
			else
			{	// Node with [0] only - empty field for this traction, no need to add force
				vfld=-1;
			}
		}
		
		// add if find a field
		if(vfld>=0)
		{	nd[nds[i]]->AddFextSpreadTask3(vfld,FTract(sign*fn[i]));
			fnorm+=fn[i];
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
	// exit if no traction law
	if(MatID()<0) return;
	
	// get tangential unit vector and length
	double dl;
	Vector t=GetTangential(&dl);
	
	// get normal and tangential COD components
	double codx=surfx[0]-surfx[1];			// above crack - below crack
	double cody=surfy[0]-surfy[1];
	double nCod=-t.y*codx+t.x*cody;			// absolute normal cod
	double tCod=t.x*codx+t.y*cody;			// absolute tangential cod
	
	// will eventually call a traction law materials
	TractionLaw *theLaw=(TractionLaw *)theMaterials[MatID()];
	theLaw->CrackTractionLaw(this,nCod,tCod,t.x,t.y,theCrack->thickness*dl);
}

// Calculate energy in the traction law for this segment if in traction
// fullEnergy true gets total energy at location of this segment or nearest traction law tip
// full Energy false gets recoverable energy at the traction law tip
// Used in J integral calculations so return energy in N/mm
double CrackSegment::TractionEnergy(Vector *crossPt,int crkTipIdx,bool fullEnergy)
{
	// if not traction law, scan toward crack tip looking for one
	if(MatID()<0)
	{	CrackSegment *closerSeg=this;
		while(TRUE)
		{	closerSeg = (crkTipIdx==START_OF_CRACK) ? closerSeg->prevSeg : closerSeg->nextSeg ;
			if(closerSeg==NULL) break;		// no traction law found near the crack tip
			
			// if find traction law, get crack tip energy (no need to interpolate)
			if(closerSeg->MatID()>=0)
				return closerSeg->SegmentTractionEnergy(fullEnergy);
		}
		
		// no traction law on this crack; should return energy due to friction if crack surfaces in contact
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
    
    // J integral (*1000 for units units J/m^2)
    if(archiver->CrackArchive(ARCH_JIntegral))
    {	*(double *)app=Jint.x*1000.;	// current crack tip J1
        app+=sizeof(double);
		if(fmobj->propagate[0])
			*(double *)app=propagationJ*1000.;		// actual energy released last time the crack grew
		else
			*(double *)app=Jint.y*1000.;		// current crack tip J2
        app+=sizeof(double);
    }
    
    // Stress Intensity Factors (/sqrt(1000) for  units MPa sqrt(m))
    if(archiver->CrackArchive(ARCH_StressIntensity))
    {	*(double *)app=sif.x*0.0316227766016838;
        app+=sizeof(double);
        *(double *)app=sif.y*0.0316227766016838;
        app+=sizeof(double);
    }
	
	// Energy Balance Results
	if(archiver->CrackArchive(ARCH_BalanceResults))
	{   *(int *)app=crackIncrements;
		app+=sizeof(double);
        *(double *)app=release*1000.;
        app+=sizeof(double);
        *(double *)app=absorb*1000.;
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

// Tell crack tip to heat itself when it propagates
void CrackSegment::StartCrackTipHeating(double growth,double tp)
{
	MaterialBase *tipMat=theMaterials[tipMatnum-1];
	double rhoH=1.0;									// should come from a material property
	double adot=tipMat->WaveSpeed(FALSE);				// in mm/secin 2D
	
	// set up rate and times (making sure proper number of steps)
	int nsteps=(int)(growth/(adot*timestep));
	if(nsteps<1) nsteps=1;
	// same as heatRate=rhoH*Jint.x*tp*adot;
	heatRate=rhoH*Jint.x*tp*growth/(timestep*(double)nsteps);	// heating rate (mm/cm)^3 W
	// discrete version of heatEndTime=mtime + growth/adot;
	heatEndTime=mtime + ((double)nsteps+.2)*timestep;			// stop heating at this time

	heating=TRUE;
}

// return heating rate if any
double CrackSegment::HeatRate(void)
{
	// return zero if off or if done
	if(!heating) return (double)0.;
	if(mtime>heatEndTime)
	{	heating=FALSE;
		return (double)0.;
	}
	
	// return the heat rate
	return heatRate;
}

// return external force (times a shape function)
Vector CrackSegment::FTract(double fni)
{	Vector fout;
	fout.x=fni*tract.x;
	fout.y=fni*tract.y;
	fout.z=0.;
	return fout;
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

	// distance crack particle to intersection place (relative to segment length)
	double t=(dxs*dxp+dys*dyp)/segLength;
	
	// if less than 0 and at first of two internal segments, return to try the next segment instead
	if(t<0. && thereIsAnotherSegement) return false;
	
	// distance surface to crack plane (relative to segment length)
	//double n=(dxs*dyp-dys*dxp)*dir/segLength;
	
	// if pretty close or negative, then do not move to plane and hope velocity fields will resolve on their own
	// or if first of two internal segments, try the other one
	//if(n<1.e-6) return !thereIsAnotherSegement;
	
	// restrict terminal segments
	if(prevSeg==NULL || nextSeg==NULL)
	{	if(t<0.) t=0.;
	}
	
	// restrict to intersect within this segment
	if(t>1.)
		t=1.;
	else if(t<-1.)
		t=-1.;
	
	// if -1 < t < 1, might want to screen out small movements using the n check which used to be above
	
	//if(x>16.)
	//	cout << "# move side " << side << " from (" << surfx[j] << "," << surfy[j] << ") to (";
	
	// move to crack plane and a little more in normal direction
	surfx[j]=x+t*dxp-1.0e-12*dir*dyp;
	surfy[j]=y+t*dyp+1.0e-12*dir*dxp;
	
	//if(x>16.)
	//	cout <<  surfx[j] << "," << surfy[j] << "), (dxp,dyp) = (" << dxp << "," << dyp << ") near (x,y) = (" << x << "," << y << ")" << endl;
	
	return true;
	
	/* old method moved along the COD See JAN-OSU-4, pg 82
	 
	 int c1=side-1;
	 int c2=1-c1;
	 double dxc=surfx[c2]-surfx[c1],dyc=surfy[c2]-surfy[c1];
	 double t;
	 
	 // separate cases depending on segment orientation
	 if(fabs(dxp)>fabs(dyp))
	 {	double slopeDiff=dxp*dyc-dyp*dxc;
	 if(fabs(slopeDiff/dxp)<1.e-6)
	 {	// close to parallel
	 double p2=dxp*dxp+dyp*dyp;
	 t=(dxp*(surfx[c1]-x)+dyp*(surfy[c1]-y))/p2;
	 }
	 else if(fabs(dxc)<1.e-8)
	 {	// vertical cod
	 t=(surfx[c1]-x)/dxp;
	 }
	 else
	 t=(dyc*(surfx[c1]-x)-dxc*(surfy[c1]-y))/slopeDiff;
	 }
	 else
	 {	double slopeDiff=dyp*dxc-dxp*dyc;
	 if(fabs(slopeDiff/dyp)<1.e-6)
	 {	// close to parallel
	 double p=dxp*dxp+dyp*dyp;
	 t=(dxp*(surfx[c1]-x)+dyp*(surfy[c1]-y))/p;
	 }
	 else if(fabs(dyc)<1.e-8)
	 {	// horizontal cod
	 t=(surfy[c1]-y)/dyp;
	 }
	 else
	 t=(dxc*(surfy[c1]-y)-dyc*(surfx[c1]-x))/slopeDiff;
	 }
	 
	 // if not in this segment, exit (unless end segment)
	 cout << " t=" << t << endl;
	 if(t<0.)
	 {	if(thereIsAnotherSegement) return false;
	 t=0.;
	 }
	 else
	 {	// keep within this segment
	 t=fmin(t,1.0);
	 }
	 
	 surfx[c1]=x+t*dxp;
	 surfy[c1]=y+t*dyp;
	 
	 return true;
	 
	 */
}

// if decide no longer need to track surfaces, can call this method to collapse surface
// to the crack plane
void CrackSegment::CollapseSurfaces(void)
{
	surfx[0]=surfx[1]=x;
	surfy[0]=surfy[1]=y;
	surfInElem[0]=surfInElem[1]=planeInElem;
}

#pragma mark ACCESSORS

// move above or below position slightly in the direction of the normal
// if it has already moved, then do not bother
Vector CrackSegment::SlightlyMoved(int side)
{
	Vector moved;
	moved.x=surfx[side-1];
	moved.y=surfy[side-1];
	moved.z=0;
	if(x!=moved.x || y!=moved.y) return moved;
	
	double dx,dy;
	if(nextSeg!=NULL)
	{	dy=x-nextSeg->x;		// -¶x
		dx=nextSeg->y-y;		// ¶y
	}
	else
	{	dy=prevSeg->x-x;		// -¶x
		dx=y-prevSeg->y;		// ¶y
	}
	
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

