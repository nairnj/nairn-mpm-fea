/********************************************************************************
    CrackHeader.cpp
    NairnMPM
    
    Created by John Nairn on Wed Apr 3 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include <fstream>

#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Exceptions/MPMTermination.hpp"
#include "Materials/MaterialBase.hpp"
#include "System/ArchiveData.hpp"
#include "Elements/ElementBase.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Cracks/ContourPoint.hpp"
#include "Cracks/CrackSegment.hpp"
#include "Cracks/CrackNode.hpp"
#include "Nodes/NodalPoint2D.hpp"
#include "MPM_Classes/MPMBase.hpp"

#ifdef HIERARCHICAL_CRACKS
#include "Cracks/CrackLeaf.hpp"
#endif

// Include to find axisymmetric Jr using Broberg method
// Comment out to use Bergkvist and Huang method
// The Broberg one appears to be much better
#define BROBERG_AS_METHOD_FOR_JR
//#define JZ_PLANAR

using namespace std; 

// class statics
double CrackHeader::codLocation;
int CrackHeader::codInterval;
double CrackHeader::bezArg[4];
double CrackHeader::bezDer[4];
int CrackHeader::warnThreeFields;
int CrackHeader::warnNodeOnCrack;
int CrackHeader::warnThreeCracks;

// globals
CrackHeader *firstCrack;		// first crack
int JGridSize = 2;				// size of J Integral contour
int JContourType = 1;			// future might try different contours
int JTerms = -1;				// number of terms in J Integral calculation (default 1 or 2 if axisymmetric)

#ifndef HIERARCHICAL_CRACKS

// extent normals (to revert to old code, set EXTENT_NORMALS to 2, and do not define HIERARCHICAL_CRACKS)
// first two must always be x=(1,0) and y=(0,1)
// rest are unnormalized P[i] = (1,enorm[i]) where enorm[i] is tangent of that angle
#if EXTENT_NORMALS == 2

#define BOUNDING_BOX_ONLY

#elif EXTENT_NORMALS == 4

// 4 Normals: x (0), y (90), +45, -45 (tangent of each angle, except for 90)
static double enorm[4] = {0.,1.,1.,-1.};

#elif EXTENT_NORMALS == 6

// 6 Normals: x (0), y (90), +60, +30, -30, -60 (tangent of each angle, except for 90)
static double enorm[6] = {0.,1.,1.732050807568877,0.5773502691896258,-0.5773502691896258,-1.732050807568877};

#endif

#endif

#pragma mark CrackHeader: Constructors and Destructors

// Constructors
CrackHeader::CrackHeader()
{	
    firstSeg=NULL;
    lastSeg=NULL;
    numberSegments=0;
	fixedCrack=FALSE;
	customContact=FALSE;
	hasTractionLaws=FALSE;
	thickness=1.0;				// for crack tip heating and tractions in mm, will default to grid thickness if set
	allowAlternate[0]=allowAlternate[1]=TRUE;
}

// Destructor
CrackHeader::~CrackHeader()
{
    CrackSegment *mseg=firstSeg,*next;
    
    while(mseg!=NULL)
    {	next=mseg->nextSeg;
        delete mseg;
        mseg=next;
    }
}

// preliminary calculations (throw CommonException on problem)
void CrackHeader::PreliminaryCrackCalcs(void)
{
    // it does not make sense unless there are two segments (and at least one line)
    if(firstSeg==NULL)
        throw CommonException("A defined crack does not have any particles","CrackHeader::PreliminaryCrackCalcs");
    else if(firstSeg->nextSeg==NULL)
        throw CommonException("All cracks must have at least two particles","CrackHeader::PreliminaryCrackCalcs");
    
	// check traction laws and set history variables
	if(hasTractionLaws)
	{	CrackSegment *scrk=firstSeg;
		while(scrk!=NULL)
		{	// check traction law
			int matid=scrk->MatID();
			if(matid>=0)
			{	if(matid>=nmat)
					throw CommonException("Crack segment with an undefined traction law material","CrackHeader::PreliminaryCrackCalcs");
				if(fixedCrack)
					throw CommonException("Fixed crack cannot have a traction law segment","CrackHeader::PreliminaryCrackCalcs");
				if(!theMaterials[matid]->isTractionLaw())
					throw CommonException("Crack segment with material that is not a traction law","CrackHeader::PreliminaryCrackCalcs");
				
				// allow traction law to have history dependent data
				scrk->SetHistoryData(theMaterials[matid]->InitHistoryData());
			}
			
			// next segment
			scrk=scrk->nextSeg;
		}
	}
		
	// check crack tip materials to be valid and to not be a traction law material
	CrackSegment *tipCrk;
	Vector crackDir;
	int crkTipIdx=START_OF_CRACK;
	while(crkTipIdx<=END_OF_CRACK)
	{	CrackTipAndDirection(crkTipIdx,&tipCrk,crackDir);
		if(tipCrk->tipMatnum>0)
		{	if(tipCrk->tipMatnum>nmat)
				throw CommonException("Crack tip material is an undefined material","CrackHeader::PreliminaryCrackCalcs");
			if(theMaterials[tipCrk->tipMatnum-1]->isTractionLaw())
				throw CommonException("Crack tip material cannot be a traction law material","CrackHeader::PreliminaryCrackCalcs");
		}
		crkTipIdx++;
	}
}

#pragma mark CrackHeader: Set up Cracks

// add new crack segment to end of the list when creating cracks
// but if same as previous location (e.g., when connecting shapes), no
// need to add this one
short CrackHeader::add(CrackSegment *cs)
{
    if(cs==NULL) return FALSE;		// not created
    if(lastSeg==NULL)
    {	firstSeg=cs;
#ifndef HIERARCHICAL_CRACKS
        CreateExtents(cs->x, cs->y);
#endif
    }
    else
	{	// no need to add a zero length segment
		if(DbleEqual(cs->x,lastSeg->x) && DbleEqual(cs->y,lastSeg->y))
		{	// but maybe want new tip material
			if(cs->tipMatnum>0)
				lastSeg->tipMatnum=cs->tipMatnum;
            if(cs->MatID()>=0)
            {   lastSeg->SetMatID(cs->MatID()+1);
                hasTractionLaws=TRUE;
            }
			return TRUE;
		}
    	lastSeg->nextSeg=cs;
		cs->prevSeg=lastSeg;
#ifndef HIERARCHICAL_CRACKS
        CheckExtents(cs->x, cs->y);
#endif
    }
    lastSeg=cs;
    numberSegments++;
	
    // determine planeInElem,surfInElem[i] for the new tip
    cs->surfInElem[0]=cs->surfInElem[1]=cs->FindElement();
	if(cs->planeInElem==0) return FALSE;
	
	// has it put traction laws on this crack
	if(cs->MatID()>=0) hasTractionLaws=TRUE;
	
   return TRUE;
}

/* add new crack segment for crack propagation 
	whichTip=0 (start) or 1 (for end)
    
    1. Add new segment at start or end
    2. Adjust firstSeg or lastSeg and nextSeg near the end points
    3. Set tipMatnum new segment to previous tip and set previous one to -1
    4. Transfer Crack tip propoerties to new crack tip
    5. Locate element for new crack tip
*/
short CrackHeader::add(CrackSegment *cs,int whichTip)
{
    CrackSegment *prevSeg;
    double dx,dy;
    int i;
    
    if(cs==NULL) return FALSE;		// not created
	
	// find the element first (problem if not in the mesh)
    cs->surfInElem[0]=cs->surfInElem[1]=cs->FindElement();
	
	// if it is not in the mesh, remove the segment and stop propagation
	if(cs->planeInElem==0)
	{	if(whichTip==END_OF_CRACK)
			lastSeg->tipMatnum=-1;
        else
			firstSeg->tipMatnum=-1;
		delete cs;
		return FALSE;
	}

    // we can assume lastSeg!=NULL because calculations will not start
    // unless all cracks have at least two particles and therefore always
    // have a lastSeg
    
    if(whichTip==END_OF_CRACK)
    {   lastSeg->nextSeg=cs;
        cs->prevSeg=lastSeg;
        prevSeg=lastSeg;
        lastSeg=cs;
    }
    else
    {   cs->nextSeg=firstSeg;
        firstSeg->prevSeg=cs;
        prevSeg=firstSeg;
        firstSeg=cs;
    }

    // transfer crack tip results to new crack tip
    cs->tipMatnum=prevSeg->tipMatnum;
    cs->Jint=prevSeg->Jint;
    cs->sif=prevSeg->sif;
    cs->propagationJ=cs->Jint.z;					// store the actual energy released
    cs->steadyState=prevSeg->steadyState;
    cs->speed=prevSeg->speed;
    for(i=0;i<3;i++)
    {   cs->potential[i]=prevSeg->potential[i];
        cs->plastic[i]=prevSeg->plastic[i];
        cs->clength[i]=prevSeg->clength[i];
    }
    cs->release=prevSeg->release;
    cs->absorb=prevSeg->absorb;
    cs->crackIncrements=prevSeg->crackIncrements;
    
    // calculate growth that formed this segment
    dx=cs->x-prevSeg->x;
    dy=cs->y-prevSeg->y;
    cs->theGrowth=sqrt(dx*dx+dy*dy);
    
    // remove original tip settings
    prevSeg->tipMatnum=-1;
    ZeroVector(&prevSeg->Jint);
    ZeroVector(&prevSeg->sif);
    
    // has it created traction law
    int tmatnum=cs->MatID();
    if(tmatnum>=0)
    {	hasTractionLaws=true;
        // history data if needed
        cs->SetHistoryData(theMaterials[tmatnum]->InitHistoryData());
    }
	
    numberSegments++;
    
#ifdef HIERARCHICAL_CRACKS
    ExtendHierarchy(cs);
#else
    CheckExtents(cs->x, cs->y);
#endif
    
   return TRUE;
}

// output crack info and evaluate contact law
void CrackHeader::Output(void)
{
	cout << "  Crack " << number << ": length = " << Length() << ", segments = " << NumberOfSegments() << ", thickness = " << thickness << " mm";
	if(fixedCrack) cout << " (fixed)";
	if(hasTractionLaws)
	{	cout << " (has traction laws)";
		// if has traction laws, must convert to frictionless
		if(!DbleEqual(crackFriction,(double)0.))
		{	customContact=TRUE;
			crackFriction=0.;
		}
	}
	cout << endl;
	if(firstSeg->tipMatnum!=-1)
	{	cout << "    start material: ";
		if(firstSeg->tipMatnum==-2)
			cout << "exterior";
		else
			cout << firstSeg->tipMatnum;
	}
	if(lastSeg->tipMatnum!=-1)
	{	cout << "    end material: ";
		if(lastSeg->tipMatnum==-2)
			cout << "exterior";
		else
			cout << lastSeg->tipMatnum;
	}
	if(firstSeg->tipMatnum!=-1 || lastSeg->tipMatnum!=-1)
		cout << endl;
	contact.CrackOutput(customContact,crackFriction,crackDn,crackDnc,crackDt,number);
	
	// save initial crack tip directions
	CrackSegment *crkTip;
	CrackTipAndDirection(START_OF_CRACK,&crkTip,initialDirection[START_OF_CRACK]);
	CrackTipAndDirection(END_OF_CRACK,&crkTip,initialDirection[END_OF_CRACK]);
	
	// future may want to read this as parameter
	SetCodLocation(1.);
	
	// check thickness
	double gridThickness=mpmgrid.GetThickness();
	if(gridThickness>0. && !DbleEqual(gridThickness,thickness))
		cout << "     +++++ WARNING: crack thickness does not match grid thickness +++++" << endl;
}

#pragma mark CrackHeader: Methods

// archive crack to file
void CrackHeader::Archive(ofstream &afile)
{
    int i=0;
    CrackSegment *mseg=firstSeg;
    
    // create space for this crack
	int recSize=archiver->GetRecordSize();
    long blen=recSize;
	char *aptr=(char *)malloc(blen);
    if(aptr==NULL)
        throw CommonException("Memory error writing crack data.","CrackHeader::Archive");
    
    // move to start and then draw
    while(mseg!=NULL)
    {	mseg->FillArchive(aptr,i);
	
 		// write this crack segment (ofstream should buffer for us)
		afile.write(aptr,blen);
		if(afile.bad())
			archiver->FileError("File error writing crack data","(results file)","CrackHeader::Archive");
			
		mseg=mseg->nextSeg;
        i++;
    }
    
    free(aptr);    
 }

// Signed area of a triangle
double CrackHeader::Triangle(double x1,double y1,double x2,double y2,double x3,double y3)
{   return(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));
}

// If contact.GetMoveOnlySurfaces() is TRUE, crack plane was already moved by
//		the surfaces; thus only need to move to midpoint and check if element has changed
// If contact.GetMoveOnlySurfaces() is FALSE, move all crack plane particles
//		using CM velocities (precalculated and stored in field[0])
// Also must recalculate extent of crack in cnear[i] and cfar[i]
short CrackHeader::MoveCrack(void)
{
	CrackSegment *scrk=firstSeg;

	// move only surfaces
	if(contact.GetMoveOnlySurfaces())
	{	while(scrk != NULL)
		{	if(!fixedCrack)
			{	// move to midpoint between upper and lower surface
				scrk->MovePosition();
				
				// did element move
				if(!scrk->FindElement()) return FALSE;
				
				// make sure surface are on correct side of the crack
				if(contact.GetPreventPlaneCrosses())
				{	if(!scrk->CheckSurfaces()) return FALSE;
				}
			}

#ifndef HIERARCHICAL_CRACKS
			// track extent
            double cx,cy;
			if(scrk==firstSeg)
			{	if(scrk->tipMatnum==EXTERIOR_CRACK)
				{   cx = 5.*scrk->x-4.*scrk->nextSeg->x;		// x1-4*(x2-x1)
					cy = 5.*scrk->y-4.*scrk->nextSeg->y;		// y1-4*(y2-y1)
                    CreateExtents(cx,cy);
				}
				else
                    CreateExtents(scrk->x,scrk->y);
			}
			else if(scrk==lastSeg && scrk->tipMatnum==EXTERIOR_CRACK)
			{	cx=5.*scrk->x-4.*scrk->prevSeg->x;		// xn+4*(xn-x(nm1))
				cy=5.*scrk->y-4.*scrk->prevSeg->y;		// yn+4*(yn-y(nm1))
                CheckExtents(cx,cy);
			}
			else
                CheckExtents(scrk->x,scrk->y);
#endif
			
			// next segments
			scrk = scrk->nextSeg;
		}
	}
	
	// move crack plane particles by CM velocity
	else
	{	int iel;
		double fn[MaxShapeNds];
		Vector cncpos;
		int j,nodeCounter;
		Vector delv,cpos,vcm;
		int nds[MaxShapeNds],numnds;
		
		// loop over crack points
		while(scrk != NULL)
		{	if(!fixedCrack)
			{	// get element and shape functinos
				iel=scrk->planeInElem-1;			// now zero based
				cpos.x=scrk->x;
				cpos.y=scrk->y;
				theElements[iel]->GetShapeFunctions(&numnds,fn,nds,&cpos,&cncpos,NULL);
				
				// initialize
				ZeroVector(&delv);
				nodeCounter=0;
				
				/*
				// renormalize shape functions in case missing some nodes
				double fnorm=0.;
				int numempty=0;
				for(j=1;j<=numnds;j++)
				{	if(nd[nds[j]]->GetCMVelocityTask8(&vcm))
					{	fnorm+=fn[j];
						AddScaledVector(&delv,&vcm,fn[j]);
					}
					else
						numempty++;
				}
				if(numempty!=0 && numempty!=numnds) ScaleVector(&delv,1./fnorm);
				*/
						
				// extrapolate to particle
				for(j=1;j<=numnds;j++)
				{	if(nd[nds[j]]->GetCMVelocityTask8(&vcm))
					{	AddScaledVector(&delv,&vcm,fn[j]);
						nodeCounter++;
					}
				}
				
				// move it or collapse it
				if(nodeCounter>0)
				{	scrk->MovePosition(timestep*delv.x,timestep*delv.y);		// in mm
					
					// did element move
					if(!scrk->FindElement()) return FALSE;
					
					// check surfaces
					if(contact.GetPreventPlaneCrosses())
					{	if(!scrk->CheckSurfaces()) return FALSE;
					}
					
					// development flag to collapse wide open cracks during cutting
					if(fmobj->dflag[0]==4 && scrk->MatID())
					{	// find COD
						double codx=scrk->surfx[0]-scrk->surfx[1];
						double cody=scrk->surfy[0]-scrk->surfy[1];
						double cod=sqrt(codx*codx+cody*cody);
						if(cod>0.75)
						{	scrk->x=(scrk->surfx[0]+scrk->surfx[1])/2.;
							scrk->y=(scrk->surfy[0]+scrk->surfy[1])/2.;
							if(!scrk->FindElement()) return FALSE;
							scrk->CollapseSurfaces();
						}
					}
				}
				//else if(scrk->MatID()<0)
				//{	// crack surface is in free space and does not have traction law
				//	scrk->CollapseSurfaces();
				//}
				else if(contact.GetPreventPlaneCrosses())
				{	// crack in free space, but check if surfaces have moved
					if(!scrk->CheckSurfaces()) return FALSE;
				}
			}


#ifndef HIERARCHICAL_CRACKS
			// track extent
            double cx,cy;
			if(scrk==firstSeg)
			{	if(scrk->tipMatnum==EXTERIOR_CRACK)
				{   cx = 5.*scrk->x-4.*scrk->nextSeg->x;		// x1-4*(x2-x1)
					cy = 5.*scrk->y-4.*scrk->nextSeg->y;		// y1-4*(y2-y1)
                    CreateExtents(cx, cy);
				}
				else
                    CreateExtents(scrk->x, scrk->y);
			}
			else if(scrk==lastSeg && scrk->tipMatnum==EXTERIOR_CRACK)
			{	cx=5.*scrk->x-4.*scrk->prevSeg->x;		// xn+4*(xn-x(nm1))
				cy=5.*scrk->y-4.*scrk->prevSeg->y;		// yn+4*(yn-y(nm1))
                CheckExtents(cx, cy);
			}
			else
                CheckExtents(scrk->x, scrk->y);
#endif
			
			// next segments
			scrk=scrk->nextSeg;
		}
	}
    
#ifdef HIERARCHICAL_CRACKS
    MoveHierarchy();
#endif
    return TRUE;
}

// Move one crack surface according to current velocities
short CrackHeader::MoveCrack(short side)
{
    CrackSegment *scrk=firstSeg;
    int iel;
    double fn[MaxShapeNds],surfaceMass;
	Vector cncpos;
    short js=side-1,nodeCounter,j;
	int numnds,nds[MaxShapeNds];
    Vector delv,cpos;
    
    // loop over crack points
    while(scrk!=NULL)
	{	if(!fixedCrack)
		{	// get element
			iel=scrk->surfInElem[js]-1;			// now zero based
			cpos.x=scrk->surfx[js];
			cpos.y=scrk->surfy[js];
			theElements[iel]->GetShapeFunctions(&numnds,fn,nds,&cpos,&cncpos,NULL);
            
			// initialize
			ZeroVector(&delv);
			surfaceMass=0;
			nodeCounter=0;
			
			// renormalize shape functions in case missing some nodes
			/*
			double fnorm=0.;
			for(j=1;j<=numnds;j++)
			{	if(nd[nds[j]]->IncrementDelvSideTask8(side,number,fn[j],&delv,&surfaceMass,scrk))
					fnorm+=fn[j];
				else
					nodeCounter++;
			}
			if(nodeCounter!=0 && nodeCounter!=numnds) ScaleVector(&delv,1./fnorm);
			
			// move it
			scrk->MoveSurfacePosition(side,timestep*delv.x,timestep*delv.y,(nodeCounter!=numnds),surfaceMass);		// in mm
			 */
			
			// extrapolate those with velocity to the particle
			for(j=1;j<=numnds;j++)
			{	if(nd[nds[j]]->IncrementDelvSideTask8(side,number,fn[j],&delv,&surfaceMass,scrk))
					nodeCounter++;
			}
			
			// delv is Sum(fi vi) = Sum(fi pi/mi) and surfaceMass = Sum(fi mi)
			if(nodeCounter>0)
			{	ScaleVector(&delv,timestep);
			}
			
			// move it (if returns true, check location of other side)
			if(scrk->MoveSurfacePosition(side,delv.x,delv.y,(nodeCounter>0),surfaceMass))		// in mm
			{	if(!scrk->FindElement(ABOVE_CRACK)) return FALSE;
			}
			
			// did surface move elements
			if(!scrk->FindElement(side)) return FALSE;
		}
            
        // on to next segement
        scrk=scrk->nextSeg;
    }
    return TRUE;
}

// Update crack tractions on any segments with traction loaws
void CrackHeader::UpdateCrackTractions(void)
{	if(!hasTractionLaws) return;
	
    // loop over crack points
	CrackSegment *scrk=firstSeg;
    while(scrk!=NULL)
	{	scrk->UpdateTractions(this);
        scrk=scrk->nextSeg;
    }
}

// called in forces step for each crack to convert crack tractions laws into external forces
void CrackHeader::TractionFext(void)
{
	if(!hasTractionLaws) return;
	CrackSegment *cs=firstSeg;
	while(cs!=NULL)
	{	cs->AddTractionFext(this);
		cs=cs->nextSeg;
	}
}

// Determine if two line-segments cross
bool CrackHeader::SegmentsCross(ContourPoint *thePt,Vector &p3,Vector &p4,Vector *crossPt)
{
    double dx,dy;
    double x1,x2,y1,y2;
    
    switch(thePt->orient)
    {	case HORIZONTAL:
            x1=thePt->node->x;
            x2=thePt->nextPoint->node->x;
            y1=thePt->node->y;
            
            // check some extents
            if(p3.y<y1 && p4.y<y1) return FALSE;
            if(p3.y>y1 && p4.y>y1) return FALSE;
            
            // find intersection
            dy=p4.y-p3.y;
            if(DbleEqual(dy,0.)) return (x1==p3.x);		// parallel lines
            dx=p4.x-p3.x;
            crossPt->y=thePt->node->y;
            crossPt->x=(crossPt->y-p3.y)*dx/dy + p3.x;
            if(x1<x2)
                return (crossPt->x>x1 && crossPt->x<=x2);
            else
                return (crossPt->x>x2 && crossPt->x<=x1);
        
        case VERTICAL:
            y1=thePt->node->y;
            y2=thePt->nextPoint->node->y;
            x1=thePt->node->x;
            
            // check some extents
            if(p3.x<x1 && p4.x<x1) return FALSE;
            if(p3.x>x1 && p4.x>x1) return FALSE;
            
            // find intersection
            dx=p4.x-p3.x;
            if(DbleEqual(dx,0.)) return (y1==p3.y);		// parallel lines
            dy=p4.y-p3.y;
            crossPt->x=thePt->node->x;
            crossPt->y=(crossPt->x-p3.x)*dy/dx + p3.y;
            if(y1<y2)
                return (crossPt->y>y1 && crossPt->y<=y2);
            else
                return (crossPt->y>y2 && crossPt->y<=y1);
        
        default:
            break;
    }
    return FALSE;
}

// load vector with initial crack tip direction
void CrackHeader::GetInitialDirection(CrackSegment *crkTip,Vector &tipDir)
{
	int whichTip=GetWhichTip(crkTip);
	tipDir.x=initialDirection[whichTip].x;
	tipDir.y=initialDirection[whichTip].y;
}

// return crack tip segment and the direction cosine vector for a crack tip
void CrackHeader::CrackTipAndDirection(int crkTipIdx,CrackSegment **crkTip,Vector &tipDir)
{
	CrackSegment *adjCrk;
	double dx=0.,dy=0.;
#ifdef _LINEAR_INTERPOLATION_
	// Find direction of last segment
	if(crkTipIdx==START_OF_CRACK)
	{	*crkTip=firstSeg;
		adjCrk=firstSeg->nextSeg;
	}
	else
	{	*crkTip=lastSeg;
		adjCrk=lastSeg->prevSeg;
	}
	dx=(*crkTip)->x-adjCrk->x;
	dy=(*crkTip)->y-adjCrk->y;
#else

#ifndef _CUBIC_INTERPOLATION_

	// Bezier curve is tangent to last segment, thus direction is same as last segment
	if(crkTipIdx==START_OF_CRACK)
	{	*crkTip=firstSeg;
		adjCrk=firstSeg->nextSeg;
	}
	else
	{	*crkTip=lastSeg;
		adjCrk=lastSeg->prevSeg;
	}
	dx=(*crkTip)->x-adjCrk->x;
	dy=(*crkTip)->y-adjCrk->y;
	
#else

	// find slope of last segment from cubic spline interpolation through the first
	//   four crack particles. Here only calculate final slope
	// formula is (1/15)*(19(S0-S1) - 5(S1-S2) + (S2-S3))
	//            (1/15 can be ignored because normalized later)
	CrackSegment *nextCrk;
	if(crkTipIdx==START_OF_CRACK)
	{	*crkTip=firstSeg;
	
		adjCrk=firstSeg;						// S0
		nextCrk=firstSeg->nextSeg;				// S1
		dx+=19.*(adjCrk->x-nextCrk->x);			// S0-S1
		dy+=19.*(adjCrk->y-nextCrk->y);
		
		adjCrk=nextCrk->nextSeg;				// S2
		if(adjCrk!=NULL)
		{	dx+=5.*(adjCrk->x-nextCrk->x);		// S2-S1
			dy+=5.*(adjCrk->y-nextCrk->y);
			
			nextCrk=adjCrk->nextSeg;			// S3
			if(nextCrk!=NULL)
			{	dx+=adjCrk->x-nextCrk->x;		// S2-S3
				dy+=adjCrk->y-nextCrk->y;
			}
		}
	}
	else
	{	*crkTip=lastSeg;
	
		adjCrk=lastSeg;						// S0
		nextCrk=lastSeg->prevSeg;				// S1
		dx+=19.*(adjCrk->x-nextCrk->x);			// S0-S1
		dy+=19.*(adjCrk->y-nextCrk->y);
		
		adjCrk=nextCrk->prevSeg;				// S2
		if(adjCrk!=NULL)
		{	dx+=5.*(adjCrk->x-nextCrk->x);		// S2-S1
			dy+=5.*(adjCrk->y-nextCrk->y);
			
			nextCrk=adjCrk->prevSeg;			// S3
			if(nextCrk!=NULL)
			{	dx+=adjCrk->x-nextCrk->x;		// S2-S3
				dy+=adjCrk->y-nextCrk->y;
			}
		}
	}
	
#endif
#endif

	// normalize
    double ds=sqrt(dx*dx+dy*dy);
    tipDir.x=dx/ds;
    tipDir.y=dy/ds;
	tipDir.z=0.;
}

// determine cod at crack tip
void CrackHeader::GetCOD(CrackSegment *crkTip,Vector &cod,bool getModes)
{
	int whichTip;
	Vector tipDir;
	
#ifdef _LINEAR_INTERPOLATION_

	CrackSegment *adjTip;
	
	if(crkTip==firstSeg)
	{	whichTip=START_OF_CRACK;
		adjTip=crkTip->nextSeg;
	}
	else
	{	whichTip=END_OF_CRACK;
		adjTip=crkTip->prevSeg;
	}
	
	// above crack - below crack, but flip if at start
	if(codLocation<=0.5)
	{	cod.x=crkTip->surfx[0]-crkTip->surfx[1];
		cod.y=crkTip->surfy[0]-crkTip->surfy[1];
	}
	else
	{	cod.x=adjTip->surfx[0]-adjTip->surfx[1];
		cod.y=adjTip->surfy[0]-adjTip->surfy[1];
	}
	if(whichTip==START_OF_CRACK)
	{	cod.x=-cod.x;
		cod.y=-cod.y;
	}
	
	// convert cod.x to shear (mode II) cod and cod.y to normal (mode I)
	if(getModes)
	{	CrackSegment *newTip;
		CrackTipAndDirection(whichTip,&newTip,tipDir);
		double codx=cod.x,cody=cod.y;
		cod.x=codx*tipDir.x+cody*tipDir.y;
		cod.y=-codx*tipDir.y+cody*tipDir.x;
	}
	
#else

	// find first four segments
	CrackSegment *seg[4];
	seg[0]=crkTip;
	if(crkTip==firstSeg)
	{	whichTip=START_OF_CRACK;
		seg[1]=crkTip->nextSeg;
		seg[2]=seg[1]->nextSeg;
		seg[3] = (seg[2]!=NULL) ? seg[2]->nextSeg : NULL;
	}
	else
	{	whichTip=END_OF_CRACK;
		seg[1]=crkTip->prevSeg;
		seg[2]=seg[1]->prevSeg;
		seg[3] = (seg[2]!=NULL) ? seg[2]->prevSeg : NULL;
	}

	// revert to linear if not enough particles
	if(seg[3]==NULL)
	{	// above crack - below crack, but flip if at start of crack
		if(codLocation<=0.5)
		{	cod.x=crkTip->surfx[0]-crkTip->surfx[1];
			cod.y=crkTip->surfy[0]-crkTip->surfy[1];
		}
		else
		{	cod.x=seg[1]->surfx[0]-seg[1]->surfx[1];
			cod.y=seg[1]->surfy[0]-seg[1]->surfy[1];
		}
		if(whichTip==START_OF_CRACK)
		{	cod.x=-cod.x;
			cod.y=-cod.y;
		}
		
		// convert cod.x to shear (mode II) cod and cod.y to normal (mode I)
		if(getModes)
		{	CrackSegment *newTip;
			CrackTipAndDirection(whichTip,&newTip,tipDir);
			double codx=cod.x,cody=cod.y;
			cod.x=codx*tipDir.x+cody*tipDir.y;
			cod.y=-codx*tipDir.y+cody*tipDir.x;
		}
		return;
	}
	
	// find raw COD
	Vector above,below;
	InterpolatePosition(ABOVE_CRACK,seg,above,false);
	InterpolatePosition(BELOW_CRACK,seg,below,false);
	cod.x=above.x-below.x;
	cod.y=above.y-below.y;
	if(whichTip==START_OF_CRACK)
	{	cod.x=-cod.x;
		cod.y=-cod.y;
	}
	
	// convert cod.x to shear (mode II) cod and cod.y to normal (mode I)
	if(getModes)
	{	//InterpolatePosition(NO_CRACK,seg,tipDir,true);
		CrackSegment *newTip;
		CrackTipAndDirection(whichTip,&newTip,tipDir);
		double codx=cod.x,cody=cod.y;
		cod.x=codx*tipDir.x+cody*tipDir.y;
		cod.y=-codx*tipDir.y+cody*tipDir.x;
	}
	
	return;
	
#endif
}

// Do cublic spline intepolations
void CrackHeader::InterpolatePosition(int surface,CrackSegment **seg,Vector &pos,bool derivative)
{
	Vector spt[4];		// interpolant points
	Vector bpt[4];		// control points
	Vector bc[4];		// bezier curve control points
	int i,j=surface-1;
	
#ifdef _CUBIC_INTERPOLATION_

	// for cubic interpolation, interpolant points are the crack poionts
	for(i=0;i<4;i++)
	{	if(surface==NO_CRACK)
		{	spt[i].x=seg[i]->x;
			spt[i].y=seg[i]->y;
		}
		else
		{	spt[i].x=seg[i]->surfx[j];
			spt[i].y=seg[i]->surfy[j];
		}
	}
	
	// calculate control points
	bpt[0].x=spt[0].x;
	bpt[0].y=spt[0].y;
	bpt[1].x=(-4.*spt[0].x+24.*spt[1].x-6.*spt[2].x+spt[3].x)/15.;
	bpt[1].y=(-4.*spt[0].y+24.*spt[1].y-6.*spt[2].y+spt[3].y)/15.;
	bpt[2].x=(spt[0].x-6.*spt[1].x+24.*spt[2].x-4.*spt[3].x)/15.;
	bpt[2].y=(spt[0].y-6.*spt[1].y+24.*spt[2].y-4.*spt[3].y)/15.;
	//bpt[3].x=spt[3].x;
	//bpt[3].y=spt[3].y;
	
#else

	// for bezier curves, crack particles are control points
	for(i=0;i<4;i++)
	{	if(surface==NO_CRACK)
		{	bpt[i].x=seg[i]->x;
			bpt[i].y=seg[i]->y;
		}
		else
		{	bpt[i].x=seg[i]->surfx[j];
			bpt[i].y=seg[i]->surfy[j];
		}
	}
	
	// calculate the interpolant points
	spt[0].x=bpt[0].x;
	spt[0].y=bpt[0].y;
	spt[1].x=(bpt[0].x+4.*bpt[1].x+bpt[2].x)/6.;
	spt[1].y=(bpt[0].y+4.*bpt[1].y+bpt[2].y)/6.;
	spt[2].x=(bpt[1].x+4.*bpt[2].x+bpt[3].x)/6.;
	spt[2].y=(bpt[1].y+4.*bpt[2].y+bpt[3].y)/6.;
	//spt[3].x=bpt[3].x;
	//spt[3].y=bpt[3].y;
	
#endif
	
	// control points for bezier curve between spt[codInterval] and spt[codInterval+1]
	// codInterval is zero ot one, thus never needs bpt[3] or spt[3]
	bc[0].x=spt[codInterval].x;
	bc[0].y=spt[codInterval].y;
	bc[1].x=(2.*bpt[codInterval].x+bpt[codInterval+1].x)/3.;
	bc[1].y=(2.*bpt[codInterval].y+bpt[codInterval+1].y)/3.;
	bc[2].x=(bpt[codInterval].x+2.*bpt[codInterval+1].x)/3.;
	bc[2].y=(bpt[codInterval].y+2.*bpt[codInterval+1].y)/3.;
	bc[3].x=spt[codInterval+1].x;
	bc[3].y=spt[codInterval+1].y;
	
	pos.x=0;
	pos.y=0.;
	if(derivative)
	{	for(i=0;i<4;i++)
		{	pos.x+=bezDer[i]*bc[i].x;
			pos.y+=bezDer[i]*bc[i].y;
		}
		double ds=sqrt(pos.x*pos.x+pos.y*pos.y);
		pos.x=-pos.x/ds;
		pos.y=-pos.y/ds;
	}
	else
	{	for(i=0;i<4;i++)
		{	pos.x+=bezArg[i]*bc[i].x;
			pos.y+=bezArg[i]*bc[i].y;
		}
	}
}

// J-integral calculation (YJG)
void CrackHeader::JIntegral(void)
{
    int gridElem,gridNode,nextNearest;
    int i,j;
    ContourPoint *crackPt,*prevPt,*nextPt;
    int crkTipIdx;
    CrackSegment *tipCrk;
    double Jx,Jy,Jx1,Jy1,Jx2,Jy2,JxAS2;
	Vector crackDir;
	bool secondTry;
	
    /* Calculate J-integrals for the ith crack tip
    */
	
	// it may try two contours at each crack tip. First try is at NearestnNode().
	// if that path crosses the crack twice, it tries a contour from the next nearest node.
	secondTry=FALSE;
	
	// each crack tip
	crkTipIdx=START_OF_CRACK;
	while(crkTipIdx<=END_OF_CRACK)
    {
		/* Task 1: find crack tip and the crack direction
        */
		CrackTipAndDirection(crkTipIdx,&tipCrk,crackDir);
        
        // done if not a crack tip
        if(tipCrk->tipMatnum<0)
		{	ZeroVector(&tipCrk->Jint);
			crkTipIdx++;
            continue;
        }
		double crackr = tipCrk->x;
		
		// block to catch problems
		NodalPoint *phantom=NULL;
		crackPt=NULL;
		try
		{
			/* Task 2: find ccw nodal points JGridSize from crack tip nodal point
				Find orientation of each line segment
			*/
			gridElem=tipCrk->planeInElem-1;
			gridNode=theElements[gridElem]->NearestNode(tipCrk->x,tipCrk->y,&nextNearest);
			if(secondTry) gridNode=nextNearest;
			int numSegs=0;
			
			// set material type based material near the tip.
			int oldnum = tipCrk->tipMatnum;
			tipCrk->FindCrackTipMaterial();
			if(tipCrk->tipMatnum != oldnum)
			{	cout << "# crack tip left material " << oldnum
						<< " and entered material " << tipCrk->tipMatnum << endl;
			}
			
			// step to the edge of the J Integral contour
			gridNode=theElements[gridElem]->NextNode(gridNode);
			for(i=1;i<JGridSize;i++)
			{   gridElem=theElements[gridElem]->Neighbor(gridNode);
				if(gridElem<0)
					throw "J integral contour at a crack tip does not fit in the grid";
				gridNode=theElements[gridElem]->NextNode(gridNode);
			}
			
			// walk around the countour (0 and 4 are half edges)
			prevPt=crackPt=new ContourPoint(nd[gridNode]);
			numSegs++;
			int numPts=JGridSize;
			double cxmin=9e99,cxmax=-9e99,cymin=9e99,cymax=-9e99;
			for(j=0;j<5;j++)
			{   gridNode=theElements[gridElem]->NextNode(gridNode);
				nextPt=new ContourPoint(nd[gridNode]);
				prevPt->SetNextPoint(nextPt);
				prevPt=nextPt;
				numSegs++;
				for(i=1;i<numPts;i++)
				{   gridElem=theElements[gridElem]->Neighbor(gridNode);
					if(gridElem<0)
						throw "J integral contour at a crack tip does not fit in the grid";
					gridNode=theElements[gridElem]->NextNode(gridNode);
					nextPt=new ContourPoint(nd[gridNode]);
					prevPt->SetNextPoint(nextPt);
					prevPt=nextPt;
					numSegs++;
				}
				
				// check corners for extent of contour rectangle
				cxmin=min(cxmin,nd[gridNode]->x);
				cxmax=max(cxmax,nd[gridNode]->x);
				cymin=min(cymin,nd[gridNode]->y);
				cymax=max(cymax,nd[gridNode]->y);

				numPts = (j==3) ? JGridSize-1 : 2*JGridSize;
			}
			// connect end to start
			prevPt->SetNextPoint(crackPt);
			
			/* Task 3: Find crack intersection with the contour
				Verify only one intersection and x-y grid
				Create phantom nodal point on the crack
			*/
			int crossCount=0;
			Vector p3,p4,crossPt,crossPt1;
			nextPt=crackPt;
			CrackSegment *startSeg=firstSeg,*endSeg;
			while(TRUE)
			{   // error if grid not alng x and y axes
				if(nextPt->orient==ANGLED)
				{	throw MPMTermination("The J Contour is not along x and y axes.",
									"CrackHeader::JIntegral");
				}
				
				//p3,p4--two end points of eack crack segment
				endSeg=firstSeg;
				p3.x=endSeg->x;
				p3.y=endSeg->y;
				endSeg=endSeg->nextSeg;
				while(endSeg!=NULL)
				{	p4.x=endSeg->x;
					p4.y=endSeg->y;
					if(SegmentsCross(nextPt,p3,p4,&crossPt1))
					{	crossCount++;
						if(crossCount>2)
						{	throw MPMTermination("More than 2 crossings between J-path and a crack",
									"CrackHeader::JIntegral");
						}
						else if(crossCount==2)
						{	if(!(DbleEqual(crossPt.x,crossPt1.x)&&DbleEqual(crossPt.y,crossPt1.y)))
							{	if(secondTry)
								{   throw MPMTermination("Two different crossings between J-path and a crack",
															"CrackHeader::JIntegral");
								}
								else
									throw "";
							}
						}
						else
						{   prevPt=nextPt;
							crossPt.x=crossPt1.x;
							crossPt.y=crossPt1.y;
							startSeg=endSeg->prevSeg;
						}
					}
					p3.x=p4.x;
					p3.y=p4.y;
					endSeg=endSeg->nextSeg;
				}
				nextPt=nextPt->nextPoint;
				if(nextPt==crackPt) break;
			}
			if(crossCount<1)
			{   throw MPMTermination("A crack does not cross its J path",
						"CrackHeader::JIntegral");
			}
			
			// find crack particle closer to the crack tipstart
			if(crkTipIdx==END_OF_CRACK) startSeg=startSeg->nextSeg;
			
			// print the contour (for debegging)
			/*
			cout << "# J Contour from node " << crackPt->node->num << ", cross at (" << crossPt.x << "," << crossPt.y << ") fraction = "
						<< prevPt->Fraction(crossPt) << " then nodes: " ;
			nextPt=prevPt->nextPoint;
			while(TRUE)
			{	cout << " " << nextPt->node->num ;
				nextPt=nextPt->nextPoint;
				if(nextPt==prevPt->nextPoint) break;
			}
			cout << endl;
			*/
			
			// insert nodal point and define start of the path
			double fract=prevPt->Fraction(crossPt);
			phantom=new NodalPoint2D(0,crossPt.x,crossPt.y);
			phantom->PrepareForFields();
			crackPt=new ContourPoint(phantom);
			crackPt->SetNextPoint(prevPt->nextPoint);
			prevPt->SetNextPoint(crackPt);
			phantom->Interpolate(prevPt->node,crackPt->nextPoint->node,fract,(tipCrk==firstSeg));				
				
			/* Task 4: Loop over all segments and evaluate J integral (from the primary term)
				Transform to crack plane and save results
			*/
			DispField *sfld1,*sfld2;
			Jx1=Jy1=0.0;			// J-integral components from the first term
			double tractionEnergy=0.,bridgingReleased=0.;
			numSegs>>=1;			// half the segments
			int dfld = (tipCrk==firstSeg) ? ABOVE_CRACK : BELOW_CRACK;		// initial field
			nextPt=crackPt;
			int count=0;			// particles in the nodal fields
			double r1 = 1.,r2 = 1.;
			while(TRUE)
			{   // J integral node1 to node2 using field dfld
				NodalPoint *node1=nextPt->node;
				NodalPoint *node2=nextPt->nextPoint->node;
				if(dfld==ABOVE_CRACK)
				{	sfld1=node1->cvf[(int)node1->above]->df;
					sfld2=node2->cvf[(int)node2->above]->df;
					count+=node2->cvf[(int)node2->above]->GetNumberPoints();
				}
				else
				{	sfld1=node1->cvf[(int)node1->below]->df;
					sfld2=node2->cvf[(int)node2->below]->df;
					count+=node2->cvf[(int)node2->below]->GetNumberPoints();
				}
				if(fmobj->IsAxisymmetric())
				{	r1 = node1->x/crackr;		// divide by a
					r2 = node2->x/crackr;
				}
				
				/* Calculate J Integral segment by segment
				   1 means the start point of the segment, 2 means the end point
				   Units work, kinetic, stress are in N/mm^2
				   Displacement gradient is dimensionless
				*/
				
				double wd1,kd1,sxx1,syy1,sxy1;
				double dudx1,dudy1,dvdx1,dvdy1;
				double wd2,kd2,sxx2,syy2,sxy2;
				double dudx2,dudy2,dvdx2,dvdy2;
				double termForJx1,termForJy1,termForJx2,termForJy2;
				double fForJx1,fForJy1,fForJx2,fForJy2;

				// segment dS and normal from ContourPoint object
				double ds=nextPt->ds;
				Vector segNorm=nextPt->norm;

				// get values of the start point
				if(sfld1!=NULL)
				{	wd1=sfld1->work;
					kd1=sfld1->kinetic;
					sxx1=sfld1->stress.xx;
					syy1=sfld1->stress.yy;
					sxy1=sfld1->stress.xy;
					dudx1=sfld1->du.x;
					dudy1=sfld1->du.y;
					dvdx1=sfld1->dv.x;
					dvdy1=sfld1->dv.y;
				}
				else
				{	wd1=0.; kd1=0.; sxx1=0.; syy1=0.; sxy1=0.;
					dudx1=0.; dudy1=0.; dvdx1=0.; dvdy1=0.;
				}

				// get values of the end point
				if(sfld2!=NULL)
				{	wd2=sfld2->work;
					kd2=sfld2->kinetic;
					sxx2=sfld2->stress.xx;
					syy2=sfld2->stress.yy;
					sxy2=sfld2->stress.xy;
					dudx2=sfld2->du.x;
					dudy2=sfld2->du.y;
					dvdx2=sfld2->dv.x;
					dvdy2=sfld2->dv.y;
				}
				else
				{	wd2=0.; kd2=0.; sxx2=0.; syy2=0.; sxy2=0.;
					dudx2=0.; dudy2=0.; dvdx2=0.; dvdy2=0.;
				}

				// calculate Jx (note that dy=segNorm.x and dx=-segNorm.y
				// or Jr is axisymmetric

				// term (ti*ui,x) (N/mm^2)
				termForJx1=(sxx1*segNorm.x+sxy1*segNorm.y)*dudx1
						  +(sxy1*segNorm.x+syy1*segNorm.y)*dvdx1;
				termForJx2=(sxx2*segNorm.x+sxy2*segNorm.y)*dudx2
						  +(sxy2*segNorm.x+syy2*segNorm.y)*dvdx2;
						  
				// [(W+K)nx-ti*ui,x] (N/mm^2)
				fForJx1=(wd1+kd1)*segNorm.x-termForJx1;
				fForJx2=(wd2+kd2)*segNorm.x-termForJx2;

				// add for two endpoints using midpoint rule
#ifdef BROBERG_AS_METHOD_FOR_JR
				Jx1+=0.5*(fForJx1 + fForJx2)*ds;		// N mm/mm^2
#else
				Jx1+=0.5*(r1*fForJx1 + r2*fForJx2)*ds;	// N mm/mm^2
#endif

				// calculate Jy (or Jz if axisymmetric)

				// term ti*ui,y
				termForJy1=(sxx1*segNorm.x+sxy1*segNorm.y)*dudy1
						  +(sxy1*segNorm.x+syy1*segNorm.y)*dvdy1;
				termForJy2=(sxx2*segNorm.x+sxy2*segNorm.y)*dudy2
						  +(sxy2*segNorm.x+syy2*segNorm.y)*dvdy2;
						  
				// [(W+K)ny-ti*ui,y]
				fForJy1=(wd1+kd1)*segNorm.y-termForJy1;
				fForJy2=(wd2+kd2)*segNorm.y-termForJy2;

				// add for two endpoints using midpoint rule
                // The r's (=r_i/a) for axisymmetric Jz integral
#ifdef JZ_PLANAR
				Jy1+=0.5*(fForJy1 + fForJy2)*ds;
#else
				Jy1+=0.5*(r1*fForJy1 + r2*fForJy2)*ds;
#endif

				// on to next segment
				numSegs--;
				if(numSegs<=0)
					dfld = (dfld==ABOVE_CRACK) ? BELOW_CRACK : ABOVE_CRACK;
				nextPt=nextPt->nextPoint;
				if(nextPt==crackPt) break;
			}
			if(count==0)
			{	nextPt=crackPt;
				NodalPoint *node1=nextPt->node;
				cout << "# J Contour intersects crack at " << node1->x << "," << node1->y << " fraction " << fract << endl;
				int dfld = (tipCrk==firstSeg) ? ABOVE_CRACK : BELOW_CRACK;		// initial field
				numSegs=2*JGridSize*4+1;
				while(TRUE)
				{   // J integral node1 to node2 using field dfld
					NodalPoint *node2=nextPt->nextPoint->node;
					if(dfld==ABOVE_CRACK)
						cout << "#  node " << node2->num << " count above " << node2->cvf[(int)node2->above]->GetNumberPoints() << endl;
					else
						cout << "#  node " << node2->num << " count below " << node2->cvf[(int)node2->below]->GetNumberPoints() << endl;
					numSegs--;
					if(numSegs<=0)
						dfld = (dfld==ABOVE_CRACK) ? BELOW_CRACK : ABOVE_CRACK;
					nextPt=nextPt->nextPoint;
					if(nextPt==crackPt) break;
				}
				throw "Section of the J Integral contour was in empty space";
			}

			/* Task 5: Evaluate J integral from the additional terms (GYJ)
				if requested */
			Jx2 = Jy2 = JxAS2 = 0.;
			if(JTerms==2)
			{   double rho,xp,yp,carea;
				double ax,ay,duxdx,duydx,duxdy,duydy;
				double vx,vy,dvxdx,dvydy,dvxdy,dvydx;
				double f2ForJx=0.,f2ForJy=0.,f2axisym=0.;
				count=0;	// number of particles within J-integral contour

				for(int p=0;p<nmpms;p++)
				{	if(theMaterials[mpm[p]->MatID()]->Rigid()) continue;
					xp=mpm[p]->pos.x;
					yp=mpm[p]->pos.y;
					if(xp>=cxmin && xp<cxmax && yp>=cymin && yp<cymax)
					{   // (xp,yp) in the contour
						count++;
						
						// Mass density g/cm^3
						rho=theMaterials[mpm[p]->MatID()]->rho;
						
						// Accelerations mm/sec^2
						Vector *acc=mpm[p]->GetAcc();
						ax=acc->x;
						ay=acc->y;
						
						// Displacement gradients (dimensionless)
						duxdx = mpm[p]->GetDuDx();
						duydy = mpm[p]->GetDvDy();
						duxdy = mpm[p]->GetDuDy();
						duydx = mpm[p]->GetDvDx();
						
						// Velocities (mm/sec)
						vx=mpm[p]->vel.x;
						vy=mpm[p]->vel.y;
						
						// Velocity gradients (1/sec)
						Tensor *velGrad=mpm[p]->GetVelGrad();
						dvxdx=velGrad->xx;
						dvydy=velGrad->yy;
						dvxdy=velGrad->xy;
						dvydx=velGrad->zz;			// yx stored in zz
						
						// increment the integrands (g/cm^3)(mm/sec^2) = N/m^3 = 1e9 N/mm^3
						f2ForJx += rho*((ax*duxdx+ay*duydx)-(vx*dvxdx+vy*dvydx)); 
						f2ForJy += rho*((ax*duxdy+ay*duydy)-(vx*dvxdy+vy*dvydy));
						
						if(fmobj->IsAxisymmetric())
						{	// in axisymmetrix z is theta direction, etheta = u/r. but w=0, az=vz=0
							// Since w=0, no change to above terms, but have some static terms for Jx=Jr only
							// Units N/(m^2 mm) = 1e6 N/mm^3
							Tensor sp = mpm[p]->ReadStressTensor();
#ifdef BROBERG_AS_METHOD_FOR_JR
							// See Broberg, Cracks and Fraction (1999), page 65
							f2axisym += rho*(sp.xx*duxdx - sp.zz*mpm[p]->GetDwDz() + sp.xy*duydx)/xp;
#else
							// Bergkvist and Huong called J3D/(a dphi)
							f2axisym += rho*(mpm[p]->GetStrainEnergy() - sp.zz*mpm[p]->GetDwDz())/crackr;
#endif
						}
					}
				}
				
				if(count==0)
					throw "J Integral contour contains no particles";
				carea=1.e-6*(cxmax-cxmin)*(cymax-cymin)/count;	// area per particle in m
				Jx2 = 1.e-3*f2ForJx*carea;				// Jx2 in N mm/mm^2 now
				Jy2 = 1.e-3*f2ForJy*carea;				// Jy2 in N mm/mm^2 now
				JxAS2 = f2axisym*carea;					// JxAS2 (for Jr in axisymmetric) in N mm/mm^2
			}
			
			/* Task 6: Subtract energy due to tractions or for cracks in contact, subtract
				energy associated with shear stress (later not yet implemented thought)
			*/
			if(hasTractionLaws)
			{	tractionEnergy=startSeg->TractionEnergy(&crossPt,crkTipIdx,true);
				bridgingReleased=startSeg->TractionEnergy(&crossPt,crkTipIdx,false);
			}
			else
			{	// set traction energy to energy due to shear if in contact (not implemented yet)
				tractionEnergy=0.;
				bridgingReleased=0.;
			}
			
			// add the two terms N mm/mm^2
			Jx = Jx1 + Jx2 - JxAS2;
			Jy = Jy1 + Jy2;

			/* Jint -- crack-axis components of dynamic J-integral
				  Jint.x is J1 in archiving and literature and is energy release rate, here
						it accounts for traction laws. Friction can be handled but not yet implemented
				  Jint.y literature J2 - needed only to convert to KI and KII (archive when propagation is off)
						physcially is J for crack growth normal to crack direction
				  Jint.z is actual energy released when the crack and traction zone propagate together
							(archived as J2 when propagation is on)
			   crackDir -- crack propagating direction cosines from above
			   Units N/mm, multiply by 1000 to get N/m = J/m^2
			*/
			tipCrk->Jint.x = Jx*crackDir.x + Jy*crackDir.y - tractionEnergy;		// Jtip or energy that will be released if crack grows
			tipCrk->Jint.y =-Jx*crackDir.y + Jy*crackDir.x;						// J2(x) - for growth normal to crack plane
			//tipCrk->Jint.y = Jx1*crackDir.x + Jy1*crackDir.y;						// J by one term (temporary)
			tipCrk->Jint.z = tipCrk->Jint.x + bridgingReleased;						// Jrel or energy released in current state
			
			// end of try block on J calculation
			secondTry=FALSE;
		}
		catch(MPMTermination term)
		{	throw term;
		}
		catch(const char *msg)
		{	// throwing "" signals to try again with the next nearest node
			if(strlen(msg)==0)
			{	secondTry=TRUE;
			}
			else
			{	cout << "# Crack No. " << number << ": "<< msg << endl;
				cout << "# J calculation and propagation will stop" << endl;
				tipCrk->tipMatnum=-1;
				ZeroVector(&tipCrk->Jint);
			}
		}
		catch( ... )
		{	throw "Unknown exception in CrackHeader::JIntegral() method";
		}
        
        /* Task 7: Release allocated objects
        */
		if(crackPt!=NULL)
		{	nextPt=crackPt;
			while(TRUE)
			{   prevPt=nextPt;
				nextPt=prevPt->nextPoint;
				delete prevPt;
				if(nextPt==crackPt || nextPt==NULL) break;
			}
		}
		if(phantom!=NULL)
			delete phantom;
		
		if(!secondTry) crkTipIdx++;
    } // end loop over crack tips
}

/* The crack should propage to (xnew,ynew)
    whichTip= 0 or 1 for start or end of crack
*/
CrackSegment *CrackHeader::Propagate(Vector &grow,int whichTip,int tractionMat)
{
    int tipMatID;
    CrackSegment *propSeg;

    // add new crack segment (new last segment)
    if(whichTip==START_OF_CRACK)
		tipMatID=firstSeg->tipMatnum;
    else
		tipMatID=lastSeg->tipMatnum;

    // create new crack segment and add to crack tip with optional traction material
    propSeg=new CrackSegment(grow.x,grow.y,tipMatID,tractionMat);
    if(!add(propSeg,whichTip))
	{	cout << "# Crack No. " << number << " left the grid; propagation was stopped" << endl;
		return NULL;
	}
	if(tractionMat>0)
	{	hasTractionLaws=TRUE;
		fmobj->SetHasTractionCracks(TRUE);
	}
	
	// return the segment
    return propSeg;
}

/* Add crack tip heating for all points in this crack
*/
void CrackHeader::CrackTipHeating(void)
{
    CrackSegment *scrk=firstSeg;
	int iel;
	double fn[MaxShapeNds];
	Vector cncpos,cpos;
	int numnds,i,nds[MaxShapeNds];
    
	// exit if no segments
    if(scrk==NULL) return;
    
    // loop over crack tips
    while(scrk!=NULL)
	{	iel=scrk->planeInElem-1;		// now zero based
		cpos.x=scrk->x;
		cpos.y=scrk->y;
		theElements[iel]->GetShapeFunctions(&numnds,fn,nds,&cpos,&cncpos,NULL);
	
        // Add crack particle heating to each node in the element
        for(i=1;i<=numnds;i++)
		{	nd[nds[i]]->fcond+=scrk->HeatRate()*fn[i];
		}
		
		// next segment
        scrk=scrk->nextSeg;
	}
}

#pragma mark GLOBAL EXTENT CRACKS

// Determine if line from particle (1) to node (2) crosses this crack
// Return ABOVE_CRACK (1), BELOW_CRACK (2), or NO_CRACK (0) and crack normal in norm
// This method uses global extents for crack, if that fails it checks all
//      segments
#ifdef HIERARCHICAL_CRACKS
short CrackHeader::FlatCrackCross(double x1,double y1,double x2,double y2,Vector *norm)
#else
short CrackHeader::CrackCross(double x1,double y1,double x2,double y2,Vector *norm)
#endif
{
    double x3,y3,x4,y4;
    CrackSegment *scrk=firstSeg;
    short cross=NO_CRACK;
    double area1,area2;
    
#ifndef HIERARCHICAL_CRACKS
    // check extents for entire crack, which may have multiple normals
    // See JANOSU-6-66
    if(fmax(x1,x2) < cnear[0]) return cross;        // Pi = (1,0)
    if(fmin(x1,x2) > cfar[0]) return cross;
    if(fmax(y1,y2) < cnear[1]) return cross;        // Pi = (0,1)
    if(fmin(y1,y2) > cfar[1]) return cross;
    
#ifndef BOUNDING_BOX_ONLY
    // remaining normals are Pi = (1,enorm[i])
    int i;
    for(i=2;i<EXTENT_NORMALS;i++)
    {   double Pia = x2-x1 + enorm[i]*(y2-y1);      // Pi.a with Pi = (1,enorm[i])
        double Pib = x1 + enorm[i]*y1;              // Pi.b
        if(Pia>0.)
        {   if(cnear[i]-Pib > Pia) return cross;
            if(cfar[i] < Pib) return cross;
        }
        else
        {   // This works for Pia=0 as well
            if(cfar[i]-Pib < Pia) return cross;
            if(cnear[i] > Pib) return cross;
        }
    }
#endif
#endif
    
    // first point
    x3=scrk->x;
    y3=scrk->y;
	if(scrk->tipMatnum==EXTERIOR_CRACK)
	{	scrk=scrk->nextSeg;
		x3-=4.*(scrk->x-x3);
		y3-=4.*(scrk->y-y3);
	}
	else
		scrk=scrk->nextSeg;

	// checking areas of various triangles
	// See JAN0048-7 for details
    while(scrk!=NULL)
    {	x4=scrk->x;
        y4=scrk->y;
		if(scrk==lastSeg)
		{	if(scrk->tipMatnum==EXTERIOR_CRACK)
			{	x4-=4.*(x3-x4);
				y4-=4.*(y3-y4);
			}
		}
        scrk=scrk->nextSeg;
        
        // check for crossing
        while(TRUE)
        {   // first two areas (123 and 124)
            area1=Triangle(x1,y1,x2,y2,x3,y3);
            area2=Triangle(x1,y1,x2,y2,x4,y4);
            
            // first area negative
            if(area1<0.)
            {	if(area2>0.)
                {   if(Triangle(x3,y3,x4,y4,x1,y1)<=0.) break;
                
                    // TRUE mean - + + (- or 0) (0 means node on crack)
                    if(Triangle(x3,y3,x4,y4,x2,y2)<=0.) goto above;
                }
                
                else if(area2==0.)
                {   if(Triangle(x3,y3,x4,y4,x1,y1)<=0.) break;
                
                    // TRUE means - 0 + 0 (node on pt 4) or - 0 + - (pt 4 between mpt and node)
                    if(Triangle(x3,y3,x4,y4,x2,y2)<=0.) goto above;
                }
            }
            
            // first area positive
            else if(area1>0.)
            {	if(area2<0.)
                {   if(Triangle(x3,y3,x4,y4,x1,y1)>=0.) break;
                    
                    // TRUE means + - - (+ or 0) (0 means node on crack)
                    if(Triangle(x3,y3,x4,y4,x2,y2)>=0.) goto below;		
                }
                
                else if(area2==0.)
                {   if(Triangle(x3,y3,x4,y4,x1,y1)>=0.) break;
                
                    // TRUE means + 0 - 0 (node on pt 4) or + 0 - + (pt 4 between mpt and node)
                    if(Triangle(x3,y3,x4,y4,x2,y2)>=0.) goto below;
                }
            }
                
            // first area zero
            else
            {	if(area2<0.)
                {   if(Triangle(x3,y3,x4,y4,x1,y1)>=0.) break;
                
                    // TRUE means 0 - - 0 (node on pt 3) or 0 - - + (pt 3 between mpt and node) 
                    if(Triangle(x3,y3,x4,y4,x2,y2)>=0.) goto below;		
                }
                
                else if(area2>0.)
                {   if(Triangle(x3,y3,x4,y4,x1,y1)<=0.) break;
                
                    // TRUE means 0 + + 0 (node on pt 3) or 0 + + - (pt 3 between mpt and node)
                    if(Triangle(x3,y3,x4,y4,x2,y2)<=0.) goto above;
                }
            }
            
            // it does not cross
            break;
			
            // toggle the setting in case there are multiple crossings
above:
            if(cross==NO_CRACK)
            {	cross=ABOVE_CRACK;
                norm->y=x3-x4;			// -¶x
                norm->x=y4-y3;			// ¶y
            }
            else
                cross=NO_CRACK;
            break;
below:
            if(cross==NO_CRACK)
            {	cross=BELOW_CRACK;
                norm->y=x3-x4;			// -¶x
                norm->x=y4-y3;			// ¶y
            }
            else
                cross=NO_CRACK;
            break;
        }
		
        // on to next segment
        x3=x4;
        y3=y4;
    }
	
    // return result
    return cross;
}

#ifndef HIERARCHICAL_CRACKS

// When adding of moving a crack, initialize the extents
// for the first segment for each normal to Pi.(cx,xy)
void CrackHeader::CreateExtents(double cx,double cy)
{
    // first normal in Pi = (1,0) so Pi.(cx,xy) = cx
    cnear[0] = cfar[0] = cx;
    
    // second normal in Pi = (0,1) so Pi.(cx,xy) = cy
    cnear[1] = cfar[1] = cy;
    
#ifndef BOUNDING_BOX_ONLY
    // remaining normals are Pi = (1,enorm[i]) so Pi.(cx,cy) = cx + enorm[i]*cy
    for(int i=2;i<EXTENT_NORMALS;i++)
        cnear[i] = cfar[i] = cx + enorm[i]*cy;
#endif

}

// When creating or movins a crack, look at new segement and
// see if need to change the extents minima or maxima
// of the term Pi.(cx,xy). cnear[i] is the global minimum
// and cfar[i] is the global maximum
void CrackHeader::CheckExtents(double cx,double cy)
{
    // first normal in Pi = (1,0) so Pi.(cx,xy) = cx
    if(cx < cnear[0]) cnear[0] = cx;
    if(cx > cfar[0]) cfar[0] = cx;
    
    // second normal in Pi = (0,1) so Pi.(cx,xy) = cy
    if(cy < cnear[1]) cnear[1] = cy;
    if(cy > cfar[1]) cfar[1] = cy;
    
#ifndef BOUNDING_BOX_ONLY
    // remaining normals are Pi = (1,enorm[i]) so Pi.(cx,cy) = cx + enorm[i]*cy
    for(int i=2;i<EXTENT_NORMALS;i++)
    {   double cij = cx + enorm[i]*cy;
        if(cij < cnear[i]) cnear[i] = cij;
        if(cij > cfar[i]) cfar[i] = cij;
    }
#endif
    
}

#endif

#pragma mark HIERARCHICAL CRACKS

#ifdef HIERARCHICAL_CRACKS

// Determine if line from particle (1) to node (2) crosses this crack
// Return ABOVE_CRACK (1), BELOW_CRACK (2), or NO_CRACK (0) and crack normal in norm
// This method uses hierarchical crack in a binary tree
short CrackHeader::CrackCross(double x1,double y1,double x2,double y2,Vector *norm)
{
    // recursive method to travese tree hierarchy
    return CrackCrossLeaf(rootLeaf,x1,y1,x2,y2,norm,NO_CRACK);
}

// Recurusive Method to process each leaf in hierarchical traversal
short CrackHeader::CrackCrossLeaf(CrackLeaf *leaf,double x1,double y1,double x2,double y2,Vector *norm,short cross)
{
    // check extents this leaf, return current cross if not in this leaf's extent
    // See JANOSU-6-66
    if(fmax(x1,x2) < leaf->cnear[0]) return cross;        // Pi = (1,0)
    if(fmin(x1,x2) > leaf->cfar[0]) return cross;
    if(fmax(y1,y2) < leaf->cnear[1]) return cross;        // Pi = (0,1)
    if(fmin(y1,y2) > leaf->cfar[1]) return cross;
    
    double Pib = x1 + y1;                       // Pi.b
    double Pia = x2 + y2 - Pib;                 // Pi.a with Pi = (1,1)
    if(Pia>0.)
    {   if(leaf->cnear[2]-Pib > Pia) return cross;
        if(leaf->cfar[2] < Pib) return cross;
    }
    else
    {   // This works for Pia=0 as well
        if(leaf->cfar[2]-Pib < Pia) return cross;
        if(leaf->cnear[2] > Pib) return cross;
    }
    
    Pib = x1 - y1;                              // Pi.b
    Pia = x2 - y2 - Pib;                        // Pi.a with Pi = (1,-1)
    if(Pia>0.)
    {   if(leaf->cnear[3]-Pib > Pia) return cross;
        if(leaf->cfar[3] < Pib) return cross;
    }
    else
    {   // This works for Pia=0 as well
        if(leaf->cfar[3]-Pib < Pia) return cross;
        if(leaf->cnear[3] > Pib) return cross;
    }
    
    // It is in extent of this leaf
    // if not terminal, go on to the child leaves
    if(!leaf->ChildrenAreSegments())
    {   CrackLeaf *child1,*child2;
        leaf->GetChildLeaves(&child1,&child2);
        cross = CrackCrossLeaf(child1,x1,y1,x2,y2,norm,cross);
        if(child2!=NULL) cross = CrackCrossLeaf(child2,x1,y1,x2,y2,norm,cross);
        return cross;
    }
    
    // Meethod 1: This code checks each segment now in a subroutine
    CrackSegment *scrk1,*scrk2;
    leaf->GetChildSegments(&scrk1,&scrk2);
    cross = CrackCrossOneSegment(scrk1,x1,y1,x2,y2,norm,cross);
    return CrackCrossOneSegment(scrk2,x1,y1,x2,y2,norm,cross);
    
    // Method 2: This code checks both segments for crossing without regard to their known extents
    // It has reached two segments. Check each one for crossing
    /*
    double x3,y3,x4,y4;
    double area1,area2;
    CrackSegment *scrk1,*scrk2;
    leaf->GetChildSegments(&scrk1,&scrk2);
    
    // first point
    CrackSegment *scrk = scrk1;
    x3=scrk1->x;
    y3=scrk1->y;
	if(scrk1==firstSeg && scrk1->tipMatnum==EXTERIOR_CRACK)
    {   x3-=4.*(scrk2->x-x3);
		y3-=4.*(scrk2->y-y3);
	}
    scrk = scrk2;
    
	// checking areas of various triangles
	// See JAN0048-7 for details
    while(scrk!=NULL)
    {	x4=scrk->x;
        y4=scrk->y;
		if(scrk==lastSeg)
		{	if(scrk->tipMatnum==EXTERIOR_CRACK)
            {	x4-=4.*(x3-x4);
                y4-=4.*(y3-y4);
            }
		}
        
        // check for crossing
        while(TRUE)
        {   // first two areas (123 and 124)
            area1=Triangle(x1,y1,x2,y2,x3,y3);
            area2=Triangle(x1,y1,x2,y2,x4,y4);
            
            // first area negative
            if(area1<0.)
            {	if(area2>0.)
                {   if(Triangle(x3,y3,x4,y4,x1,y1)<=0.) break;
                    
                    // TRUE mean - + + (- or 0) (0 means node on crack)
                    if(Triangle(x3,y3,x4,y4,x2,y2)<=0.) goto above;
                }
                    
                else if(area2==0.)
                {   if(Triangle(x3,y3,x4,y4,x1,y1)<=0.) break;
                    
                    // TRUE means - 0 + 0 (node on pt 4) or - 0 + - (pt 4 between mpt and node)
                    if(Triangle(x3,y3,x4,y4,x2,y2)<=0.) goto above;
                }
            }
            
            // first area positive
            else if(area1>0.)
            {	if(area2<0.)
                {   if(Triangle(x3,y3,x4,y4,x1,y1)>=0.) break;
                    
                    // TRUE means + - - (+ or 0) (0 means node on crack)
                    if(Triangle(x3,y3,x4,y4,x2,y2)>=0.) goto below;		
                }
                    
                else if(area2==0.)
                {   if(Triangle(x3,y3,x4,y4,x1,y1)>=0.) break;
                    
                    // TRUE means + 0 - 0 (node on pt 4) or + 0 - + (pt 4 between mpt and node)
                    if(Triangle(x3,y3,x4,y4,x2,y2)>=0.) goto below;
                }
            }
            
            // first area zero
            else
            {	if(area2<0.)
                {   if(Triangle(x3,y3,x4,y4,x1,y1)>=0.) break;
                    
                    // TRUE means 0 - - 0 (node on pt 3) or 0 - - + (pt 3 between mpt and node) 
                    if(Triangle(x3,y3,x4,y4,x2,y2)>=0.) goto below;		
                }
                    
                else if(area2>0.)
                {   if(Triangle(x3,y3,x4,y4,x1,y1)<=0.) break;
                    
                    // TRUE means 0 + + 0 (node on pt 3) or 0 + + - (pt 3 between mpt and node)
                    if(Triangle(x3,y3,x4,y4,x2,y2)<=0.) goto above;
                }
            }
            
            // it does not cross
            break;
			
            // toggle the setting in case there are multiple crossings
        above:
            if(cross==NO_CRACK)
            {	cross=ABOVE_CRACK;
                norm->y=x3-x4;			// -¶x
                norm->x=y4-y3;			// ¶y
            }
            else
                cross=NO_CRACK;
            break;
        below:
            if(cross==NO_CRACK)
            {	cross=BELOW_CRACK;
                norm->y=x3-x4;			// -¶x
                norm->x=y4-y3;			// ¶y
            }
            else
                cross=NO_CRACK;
            break;
        }
		
        // on to next segment
        if(scrk!=scrk2) break;
        scrk=scrk->nextSeg;
        x3=x4;
        y3=y4;
    }
    
    return cross;
     */
}

// Deepest Method to process one line segment in crack from particle scrk1 to scrk1->nextSeg
// Checks extents of that segment first. If fails, finally do line segment crossing algorithm
short CrackHeader::CrackCrossOneSegment(CrackSegment *scrk1,double x1,double y1,double x2,double y2,Vector *norm,short cross)
{
    // next segment, exit if none (i.e., terminal particle)
    CrackSegment *scrk2 = scrk1->nextSeg;
    if(scrk2 == NULL) return cross;
    
    // check this segment's extents
    // See JANOSU-6-66
    if(fmax(x1,x2) < scrk1->cnear[0]) return cross;        // Pi = (1,0)
    if(fmin(x1,x2) > scrk1->cfar[0]) return cross;
    if(fmax(y1,y2) < scrk1->cnear[1]) return cross;        // Pi = (0,1)
    if(fmin(y1,y2) > scrk1->cfar[1]) return cross;
    
    double Pib = x1 + y1;                       // Pi.b
    double Pia = x2 + y2 - Pib;                 // Pi.a with Pi = (1,1)
    if(Pia>0.)
    {   if(scrk1->cnear[2]-Pib > Pia) return cross;
        if(scrk1->cfar[2] < Pib) return cross;
    }
    else
    {   // This works for Pia=0 as well
        if(scrk1->cfar[2]-Pib < Pia) return cross;
        if(scrk1->cnear[2] > Pib) return cross;
    }
    
    Pib = x1 - y1;                              // Pi.b
    Pia = x2 - y2 - Pib;                        // Pi.a with Pi = (1,-1)
    if(Pia>0.)
    {   if(scrk1->cnear[3]-Pib > Pia) return cross;
        if(scrk1->cfar[3] < Pib) return cross;
    }
    else
    {   // This works for Pia=0 as well
        if(scrk1->cfar[3]-Pib < Pia) return cross;
        if(scrk1->cnear[3] > Pib) return cross;
    }
    
    // No must check for crossing
    double x3,y3,x4,y4;
    double area1,area2;
    
    // first point
    x3=scrk1->x;
    y3=scrk1->y;
	if(scrk1==firstSeg && scrk1->tipMatnum==EXTERIOR_CRACK)
    {   x3-=4.*(scrk2->x-x3);
		y3-=4.*(scrk2->y-y3);
	}
    
    // second point
    x4=scrk2->x;
    y4=scrk2->y;
    if(scrk2==lastSeg)
    {	if(scrk2->tipMatnum==EXTERIOR_CRACK)
        {	x4-=4.*(x3-x4);
            y4-=4.*(y3-y4);
        }
    }
        
	// checking areas of various triangles
	// See JAN0048-7 for details
    while(TRUE)
    {   // first two areas (123 and 124)
        area1=Triangle(x1,y1,x2,y2,x3,y3);
        area2=Triangle(x1,y1,x2,y2,x4,y4);
        
        // first area negative
        if(area1<0.)
        {	if(area2>0.)
            {   if(Triangle(x3,y3,x4,y4,x1,y1)<=0.) break;
            
                // TRUE mean - + + (- or 0) (0 means node on crack)
                if(Triangle(x3,y3,x4,y4,x2,y2)<=0.) goto above;
            }
            
            else if(area2==0.)
            {   if(Triangle(x3,y3,x4,y4,x1,y1)<=0.) break;
            
                // TRUE means - 0 + 0 (node on pt 4) or - 0 + - (pt 4 between mpt and node)
                if(Triangle(x3,y3,x4,y4,x2,y2)<=0.) goto above;
            }
        }
        
        // first area positive
        else if(area1>0.)
        {	if(area2<0.)
            {   if(Triangle(x3,y3,x4,y4,x1,y1)>=0.) break;
                
                // TRUE means + - - (+ or 0) (0 means node on crack)
                if(Triangle(x3,y3,x4,y4,x2,y2)>=0.) goto below;		
            }
                
            else if(area2==0.)
            {   if(Triangle(x3,y3,x4,y4,x1,y1)>=0.) break;
                
                // TRUE means + 0 - 0 (node on pt 4) or + 0 - + (pt 4 between mpt and node)
                if(Triangle(x3,y3,x4,y4,x2,y2)>=0.) goto below;
            }
        }
        
        // first area zero
        else
        {	if(area2<0.)
            {   if(Triangle(x3,y3,x4,y4,x1,y1)>=0.) break;
                
                // TRUE means 0 - - 0 (node on pt 3) or 0 - - + (pt 3 between mpt and node) 
                if(Triangle(x3,y3,x4,y4,x2,y2)>=0.) goto below;		
            }
                
            else if(area2>0.)
            {   if(Triangle(x3,y3,x4,y4,x1,y1)<=0.) break;
                
                // TRUE means 0 + + 0 (node on pt 3) or 0 + + - (pt 3 between mpt and node)
                if(Triangle(x3,y3,x4,y4,x2,y2)<=0.) goto above;
            }
        }
        
        // it does not cross
        break;
        
        // toggle the setting in case there are multiple crossings
    above:
        if(cross==NO_CRACK)
        {	cross=ABOVE_CRACK;
            norm->y=x3-x4;			// -¶x
            norm->x=y4-y3;			// ¶y
        }
        else
            cross=NO_CRACK;
        break;
    below:
        if(cross==NO_CRACK)
        {	cross=BELOW_CRACK;
            norm->y=x3-x4;			// -¶x
            norm->x=y4-y3;			// ¶y
        }
        else
            cross=NO_CRACK;
        break;
    }
		    
    return cross;
}

// Compare Hierarchical crack crossing result to flat one
// A call can be inserted after each call to CrackCrossing to see if hierarchical method gets a different
//   answer than the flat one
// The option at the end revert to flat results
void CrackHeader::CFFlatCrossing(double x1,double y1,double x2,double y2,Vector *norm,short *vfld,int pnum,int nodenum)
{
    Vector flatNorm;
    short svfld = *vfld;
    short flatVfld = FlatCrackCross(x1,y1,x2,y2,&flatNorm);
    if(flatVfld!=svfld || (svfld!=NO_CRACK && (!DbleEqual(norm->x,flatNorm.x) || !DbleEqual(norm->y,flatNorm.y))) )
    {   cout << "# Discrepancy for cross from partcle ";
        if(pnum>0)
            cout << pnum;
        else
            cout << "at (" << x1 << "," << y1 << ")";
        cout << " to node ";
        if(nodenum>0)
            cout << nodenum;
        else
            cout << "at (" << x2 << "," << y2 << ")";
        cout << endl;
        cout << "#   Hier: " << *vfld << ", " << norm->x << ", " << norm->y << endl;
        cout << "#   Flat: " << flatVfld << ", " << flatNorm.x << ", " << flatNorm.y << endl;
        
        // option to revert to flat calculation
        //*vfld = flatVfld;
        //norm->x = flatNorm.x;
        //norm->y = flatNorm.y;
    }
}

// When crack is first created at start of calculations, create all the CrackLeaf
// objects needed to describe the crack as a binary tree starting from
// the rootleaf
void CrackHeader::CreateHierarchy(void)
{
    CrackSegment *scrk1 = firstSeg;
    CrackSegment *scrk2 = scrk1->nextSeg;           // firstSeg always exists
    CrackLeaf *firstLeaf=NULL,*lastLeaf,*leaf;
    int numLeaves=0;
    
    // consider line segments in pairs (or could be only 1 if scrk2->nextSeg==NULL)
    //  scrk1 to scrk2 and scrk2 to scrk2->nextSeg
    while(scrk2!=NULL)
    {   // set extents of each segment
        scrk1->CreateSegmentExtents(scrk1==firstSeg);
        scrk2->CreateSegmentExtents(FALSE);
        
        // make a leaf
        leaf = new CrackLeaf(scrk1,scrk2);
        numLeaves++;
        if(firstLeaf==NULL)
            firstLeaf = leaf;
        else
            lastLeaf->nextLeaf = leaf;
        
        // on to next one
        scrk1 = scrk2->nextSeg;
        if(scrk1==NULL) break;
        scrk2 = scrk1->nextSeg;
        lastLeaf = leaf;
    }
	
    // trace up the tree as long as more than 1 leaf is present
    while(numLeaves>1)
    {   CrackLeaf *leaf1 = firstLeaf;
        CrackLeaf *leaf2 = leaf1->nextLeaf;
        firstLeaf = NULL;
        numLeaves = 0;
        
        while(leaf1!=NULL)
        {   // make parent leaf
            leaf = new CrackLeaf(leaf1,leaf2);
            numLeaves++;
            if(firstLeaf==NULL)
                firstLeaf = leaf;
            else
                lastLeaf->nextLeaf = leaf;
            
            // on to next one
            if(leaf2==NULL) break;
            leaf1 = leaf2->nextLeaf;
            if(leaf1!=NULL)
            {   leaf2 = leaf1->nextLeaf;
                lastLeaf = leaf;
            }
        }
    }
    
    // all done, but root is the most recent firstLeaf
    rootLeaf = firstLeaf;
}

// When crack is first created at start of calculations, create all the CrackLeaf
// objects needed to describe the crack as a binary tree starting from
// the rootleaf
void CrackHeader::MoveHierarchy(void)
{	
    CrackSegment *scrk1 = firstSeg;
    CrackSegment *scrk2 = scrk1->nextSeg;           // firstSeg always exists
    
    // consider line segments in pairs (or could be only 1 if scrk2->nextSeg==NULL)
    //  scrk1 to scrk2 and scrk2 to scrk2->nextSeg
    while(scrk2!=NULL)
    {   // set extents of each segment
        scrk1->CreateSegmentExtents(scrk1==firstSeg);
        scrk2->CreateSegmentExtents(FALSE);
        
        // on to next one
        scrk1 = scrk2->nextSeg;
        if(scrk1==NULL) break;
        scrk2 = scrk1->nextSeg;
    }
    
    // trace up the tree to the root leaf
    CrackLeaf *firstLeaf = firstSeg->parent;
    while(firstLeaf!=NULL)
    {   CrackLeaf *leaf = firstLeaf;
        
        // recalculate extents in this level
        while(leaf!=NULL)
        {   leaf->GetLeafExtents();
            leaf = leaf->nextLeaf;
        }
        
        // on to parent level
        firstLeaf = firstLeaf->parent;
    }
}

// When crack propagates, a new segement (here cs) is added at the start (if cs==firstSeg) or at
// then end (if cs==lastSeg). The crack hierarchy needs to accomodate the new crack segment
void CrackHeader::ExtendHierarchy(CrackSegment *cs)
{
    // add segment at the beginning of the crack
    if(cs == firstSeg)
    {   // shift leaves
        CrackSegment *scrk1 = firstSeg;
        CrackSegment *scrk2 = scrk1->nextSeg;
        
        while(scrk2!=NULL)
        {   // shift leaf one spot to the left or create new one at the end
            if(scrk2->nextSeg!=NULL)
            {   CrackLeaf *leaf = scrk2->parent;
                leaf->ShiftLeft(scrk1);
            }
            else
            {   // needs new leaf on the end, when new crack has even number of particles
                ExtendHierarchy(scrk2);
                return;
            }
        
            // on to next one
            scrk1 = scrk2->nextSeg;
            if(scrk1==NULL) break;
            scrk2 = scrk1->nextSeg;
        }
    }
    
    // add a segment at the end of the crack
    else
    {   // the current last leaf is parent of penultimate crack particle
        CrackSegment *oldLast = cs->prevSeg;
        CrackSegment *oldPenSeg = oldLast->prevSeg;
        CrackLeaf *leaf = oldPenSeg->parent;
        
        // try to add to previous last leaf, but when new crack has even number of particles
        // it will require a new leaf
        CrackLeaf *newLeaf = leaf->AddSegment(oldLast);
        
        // keeping adding leaves until find one that has room or reach the end
        while(newLeaf !=NULL)
        {   leaf->nextLeaf = newLeaf;
            if(leaf == rootLeaf)
            {   rootLeaf = new CrackLeaf(leaf,newLeaf);
                break;
            }
            leaf = leaf->parent;
            newLeaf = leaf->AddLeaf(newLeaf);
        }
    }
    
    // Once all leaves rearrange, reuse Move code to update extents
    MoveHierarchy();
    //Describe();
    //rootLeaf->DescribeSegments(0);
}

#endif

#pragma mark ACCESSORS

// number of full cracks this crack to end of list
// so far only used when printing out problem
int CrackHeader::Count(void)
{
    int ncracks=1;
    CrackHeader *next=(CrackHeader *)nextObject;
    
    while(next!=NULL)
    {	ncracks++;
        next=(CrackHeader *)next->GetNextObject();
    }
    return ncracks;
}

// Find length of current crack
double CrackHeader::Length(void)
{
    CrackSegment *scrk=firstSeg;
    double length=0.,lastx,lasty,dist2,x,y;
    
    // first point
    if(scrk==NULL) return length;
    lastx=(scrk->surfx[0]+scrk->surfx[1])/2.;
    lasty=(scrk->surfy[0]+scrk->surfy[1])/2.;
    scrk=scrk->nextSeg;
    
    // loop over crack points
    while(scrk!=NULL)
    {	// add segment length
    	x=(scrk->surfx[0]+scrk->surfx[1])/2.;
        y=(scrk->surfy[0]+scrk->surfy[1])/2.;
        dist2=(x-lastx)*(x-lastx) + (y-lasty)*(y-lasty);
        length+=sqrt(dist2);
        lastx=x;
        lasty=y;
        
        // next segment
        scrk=scrk->nextSeg;
    }
    return length;
}

// Count Segments
int CrackHeader::NumberOfSegments(void) { return numberSegments; }

// crack number
void CrackHeader::SetNumber(int newNum) { number=newNum; }
int CrackHeader::GetNumber(void) { return number; }

// make it a fixed crack
void CrackHeader::SetFixedCrack(int fixSetting) { fixedCrack=fixSetting; }

// set custom friction coefficient and interface parameters
void CrackHeader::SetFriction(double frict)
{	customContact=TRUE;
	crackFriction=frict;
}
void CrackHeader::SetDn(double Dn)
{	customContact=TRUE;
	crackDn=Dn;
	if(crackFriction<10.) crackFriction=11.;
}
void CrackHeader::SetDnc(double Dnc)
{	customContact=TRUE;
	crackDnc=Dnc;
	if(crackFriction<10.) crackFriction=11.;
}
void CrackHeader::SetDt(double Dt)
{	customContact=TRUE;
	crackDt=Dt;
	if(crackFriction<10.) crackFriction=11.;
}

// set default values
void CrackHeader::SetContact(double frict,double Dn,double Dnc,double Dt)
{	customContact=FALSE;
	crackFriction=frict;
	crackDn=Dn;
	crackDnc=Dnc;
	crackDt=Dt;
}

// given whichTip find crack tip segment
CrackSegment *CrackHeader::GetCrackTip(int whichTip)
{	return (whichTip==START_OF_CRACK) ? firstSeg : lastSeg ;
}

// given whichTip find crack tip segment
CrackSegment *CrackHeader::GetAdjToCrackTip(int whichTip)
{	return (whichTip==START_OF_CRACK) ? firstSeg->nextSeg : lastSeg->prevSeg ;
}

// given crack tip segment find which tip
int CrackHeader::GetWhichTip(CrackSegment *crkTip)
{	return (crkTip==firstSeg) ? START_OF_CRACK  : END_OF_CRACK ;
}

// return if any segment has traction law material
bool CrackHeader::GetHasTractionLaws(void) { return hasTractionLaws; }

// Deterimine what needs this crack tip has to do propagation calculation
// Return 0 (no need), NEED_JANDK, or NEED_J
int CrackHeader::CriterionNeeds(void)
{
	int thisCrackNeeds=0;
	int crkTipIdx;
	CrackSegment *tipCrk;
	
	// check crack tips with propagation possible
    for(crkTipIdx=START_OF_CRACK;crkTipIdx<=END_OF_CRACK;crkTipIdx++)
    {	// get tip info
		tipCrk=GetCrackTip(crkTipIdx);
		
		// if not tip material, then not propagating there
		if(tipCrk->tipMatnum<0) continue;
		
		// check crack tip material criterion (and alternate criterion if there)
		thisCrackNeeds|=theMaterials[tipCrk->tipMatnum-1]->CriterionNeeds(0);
		if(GetAllowAlternate(crkTipIdx))
			thisCrackNeeds|=theMaterials[tipCrk->tipMatnum-1]->CriterionNeeds(1);
	}
	
	// return the result
	return thisCrackNeeds;
}

// load vector with initial crack tip direction
bool CrackHeader::GetAllowAlternate(int crkTipIdx) { return allowAlternate[crkTipIdx]; }
void CrackHeader::SetAllowAlternate(int crkTipIdx,bool setting) { allowAlternate[crkTipIdx]=setting; }

// describe velocity field
void CrackHeader::Describe(void)
{	cout << "# crack#=" << number << " thickness=" << thickness << " custom contact=" << customContact << " has traction=" << hasTractionLaws << endl;
    cout << "#    from (" << firstSeg->x << "," << firstSeg->y << ") to (" << lastSeg->x << "," << lastSeg->y << "), length = " << Length() ;
    cout << ", segments = " << NumberOfSegments() << endl;
}

// crack thickness for planar problems
void CrackHeader::SetThickness(double thick) { thickness = thick; }
double CrackHeader::GetThickness(void) { return thickness; }
double *CrackHeader::GetThicknessPtr(void) { return &thickness; }

#pragma mark CrackHeader: Class methods

#include "System/ArchiveData.hpp"
/*
	Check for contact on crack surfaces.
	
	As explained in nodalVelBC::GridMomentumConditions(), the
	first pass has makeCopy=TRUE which copies initial
	velocites/momenta. On second pass, the acceleration/force
	is adjusted to be consistent with the nodal velocities
*/
void CrackHeader::ContactConditions(int makeCopy)
{
	if(firstCrack==NULL) return;
	
	// look for crack contact in these momenta
	int i;
	if(makeCopy)
	{	for(i=1;i<=nnodes;i++)
		{	nd[i]->CrackContact(makeCopy,FALSE,0.);
		}
	}
	else
		CrackNode::ContactOnKnownNodes();
}

// Find location for spline interpolation on crack surfaces
void CrackHeader::SetCodLocation(double t)
{
	if(t<0.) t=0.;
	if(t>2.) t=2.;
	
	if(t<=1.)
		codInterval=0;
	else
	{	codInterval=1;
		t-=1.;
	}
	codLocation=t;
	
	bezArg[0]=(1.-t)*(1.-t)*(1.-t);
	bezArg[1]=3.*t*(1.-t)*(1.-t);
	bezArg[2]=3.*t*t*(1.-t);
	bezArg[3]=t*t*t;
	
	bezDer[0]=-3.*(1.-t)*(1.-t);
	bezDer[1]=3.*(1.-t)*(1.-t) - 6.*t*(1.-t);
	bezDer[2]=6.*t*(1-t) - 3.*t*t;
	bezDer[3]=3.*t*t;
}
