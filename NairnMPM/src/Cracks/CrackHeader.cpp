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
CrackHeader *firstCrack;	// first crack
int JGridSize=2;			// size of J Integral contour
int JContourType=1;			// future might try different contours
int JTerms=1;				// number of terms in J Integral calculation

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
				scrk->SetHistoryData(theMaterials[matid]->MaterialData());
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
        xmin=xmax=cs->x;
        ymin=ymax=cs->y;
    }
    else
	{	// no need to add a zero length segment
		if(DbleEqual(cs->x,lastSeg->x) && DbleEqual(cs->y,lastSeg->y))
		{	// but maybe want new tip material
			if(cs->tipMatnum>0)
				lastSeg->tipMatnum=cs->tipMatnum;
			return TRUE;
		}
    	lastSeg->nextSeg=cs;
		cs->prevSeg=lastSeg;
        xmin=fmin(xmin,cs->x);
        xmax=fmax(xmax,cs->x);
        ymin=fmin(ymin,cs->y);
        ymax=fmax(ymax,cs->y);
    }
    lastSeg=cs;
    numberSegments++;
	
    // determine planeInElem,surfInElem[i] for the new tip
    cs->surfInElem[0]=cs->surfInElem[1]=cs->FindElement();
	if(cs->planeInElem==0) return FALSE;
	
	// has it put traction laws on this crack
	if(cs->MatID()>=0) hasTractionLaws=true;
	
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

    if(lastSeg==NULL)
    {	// should never happen
    	lastSeg=firstSeg=cs;
        xmin=xmax=cs->x;
        ymin=ymax=cs->y;
    }
    else
    {	if(whichTip==END_OF_CRACK)
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
        xmin=fmin(xmin,cs->x);
        xmax=fmax(xmax,cs->x);
        ymin=fmin(ymin,cs->y);
        ymax=fmax(ymax,cs->y);
        
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
			cs->SetHistoryData(theMaterials[tmatnum]->MaterialData());
		}
	
    }
    numberSegments++;
    
    return TRUE;
}

// output crack info and evaluate contact law
void CrackHeader::Output(void)
{
	cout << "  Crack " << number << ": length = " << Length() << ", thickness = " << thickness << " mm";
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
short CrackHeader::MoveCrack(void)
{
	CrackSegment *scrk=firstSeg;

	// move only surfaces
	if(contact.GetMoveOnlySurfaces())
	{	while(scrk!=NULL)
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
			
			// track extent
			if(scrk==firstSeg)
			{	if(scrk->tipMatnum==EXTERIOR_CRACK)
				{   xmin=xmax=5.*scrk->x-4.*scrk->nextSeg->x;		// x1-4*(x2-x1)
					ymin=ymax=5.*scrk->y-4.*scrk->nextSeg->y;		// y1-4*(y2-y1)
				}
				else
				{   xmin=xmax=scrk->x;
					ymin=ymax=scrk->y;
				}
			}
			else if(scrk==lastSeg && scrk->tipMatnum==EXTERIOR_CRACK)
			{	double cx=5.*scrk->x-4.*scrk->prevSeg->x;		// xn+4*(xn-x(nm1))
				double cy=5.*scrk->y-4.*scrk->prevSeg->y;		// yn+4*(yn-y(nm1))
			    xmin=fmin(xmin,cx);
				xmax=fmax(xmax,cx);
				ymin=fmin(ymin,cy);
				ymax=fmax(ymax,cy);
			}
			else
			{   xmin=fmin(xmin,scrk->x);
				xmax=fmax(xmax,scrk->x);
				ymin=fmin(ymin,scrk->y);
				ymax=fmax(ymax,scrk->y);
			}
			
			// next segments
			scrk=scrk->nextSeg;
		}
	}
	
	// move crack plane particles by CM velocity
	else
	{	long iel;
		double fn[MaxShapeNds];
		Vector cncpos;
		int j,nodeCounter;
		Vector delv,cpos,vcm;
		int nds[MaxShapeNds],numnds;
		
		// loop over crack points
		while(scrk!=NULL)
		{	if(!fixedCrack)
			{	// get element and shape functinos
				iel=scrk->planeInElem-1;			// now zero based
				cpos.x=scrk->x;
				cpos.y=scrk->y;
				theElements[iel]->GetShapeFunctions(&numnds,fn,nds,&cpos,&cncpos);
				
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
			
			// track extent
			if(scrk==firstSeg)
			{	if(scrk->tipMatnum==EXTERIOR_CRACK)
				{   xmin=xmax=5.*scrk->x-4.*scrk->nextSeg->x;		// x1-4*(x2-x1)
					ymin=ymax=5.*scrk->y-4.*scrk->nextSeg->y;		// y1-4*(y2-y1)
				}
				else
				{   xmin=xmax=scrk->x;
					ymin=ymax=scrk->y;
				}
			}
			else if(scrk==lastSeg && scrk->tipMatnum==EXTERIOR_CRACK)
			{	double cx=5.*scrk->x-4.*scrk->prevSeg->x;		// xn+4*(xn-x(nm1))
				double cy=5.*scrk->y-4.*scrk->prevSeg->y;		// yn+4*(yn-y(nm1))
			    xmin=fmin(xmin,cx);
				xmax=fmax(xmax,cx);
				ymin=fmin(ymin,cy);
				ymax=fmax(ymax,cy);
			}
			else
			{   xmin=fmin(xmin,scrk->x);
				xmax=fmax(xmax,scrk->x);
				ymin=fmin(ymin,scrk->y);
				ymax=fmax(ymax,scrk->y);
			}
			
			// next segments
			scrk=scrk->nextSeg;
		}
	}
	
    return TRUE;
}

// Move one crack surface according to current velocities
short CrackHeader::MoveCrack(short side)
{
    CrackSegment *scrk=firstSeg;
    long iel;
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
			theElements[iel]->GetShapeFunctions(&numnds,fn,nds,&cpos,&cncpos);
            
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
    long gridElem,gridNode,nextNearest;
    long i,j;
    ContourPoint *crackPt,*prevPt,*nextPt;
    int crkTipIdx;
    CrackSegment *tipCrk;
    double Jx,Jy,Jx1,Jy1,Jx2,Jy2;
	Vector crackDir;
	bool secondTry;
	
    /* Calculate J-integrals for the ith crack tip
    */
	
	// it may try to contours at each crack tip. First try is at NearestnNode().
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
				cymax=max(cxmax,nd[gridNode]->y);

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
			phantom=new NodalPoint2D((long)0,crossPt.x,crossPt.y);
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

				// calculate Jx

				// term (ti*ui,x) (N/mm^2)
				termForJx1=(sxx1*segNorm.x+sxy1*segNorm.y)*dudx1
						  +(sxy1*segNorm.x+syy1*segNorm.y)*dvdx1;
				termForJx2=(sxx2*segNorm.x+sxy2*segNorm.y)*dudx2
						  +(sxy2*segNorm.x+syy2*segNorm.y)*dvdx2;
						  
				// [(W+K)nx-ti*ui,x] (N/mm^2)
				fForJx1=(wd1+kd1)*segNorm.x-termForJx1;
				fForJx2=(wd2+kd2)*segNorm.x-termForJx2;

				// add for two endpoints using midpoint rule
				Jx1+=0.5*(fForJx1+fForJx2)*ds;	// N mm/mm^2

				// calculate Jy

				// term ti*ui,y
				termForJy1=(sxx1*segNorm.x+sxy1*segNorm.y)*dudy1
						  +(sxy1*segNorm.x+syy1*segNorm.y)*dvdy1;
				termForJy2=(sxx2*segNorm.x+sxy2*segNorm.y)*dudy2
						  +(sxy2*segNorm.x+syy2*segNorm.y)*dvdy2;
						  
				// [(W+K)ny-ti*ui,y]
				fForJy1=(wd1+kd1)*segNorm.y-termForJy1;
				fForJy2=(wd2+kd2)*segNorm.y-termForJy2;

				// add for two endpoints using midpoint rule
				Jy1+=0.5*(fForJy1+fForJy2)*ds; 

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
			Jx2=Jy2=0.;
			if(JTerms==2)
			{   double rho,xp,yp,carea;
				double ax,ay,duxdx,duydx,duxdy,duydy;
				double vx,vy,dvxdx,dvydy,dvxdy,dvydx;
				double f2ForJx=0.,f2ForJy=0.;
				count=0;	// number of particles within J-integral contour

				for(long p=0;p<nmpms;p++)
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
						Tensor *ep=mpm[p]->GetStrainTensor();
						duxdx=ep->xx;
						duydy=ep->yy;
						duxdy=mpm[p]->GetDuDy();
						duydx=mpm[p]->GetDvDx();
						// Velocities (mm/sec)
						vx=mpm[p]->vel.x;
						vy=mpm[p]->vel.y;
						// Velocity gradients
						Tensor *velGrad=mpm[p]->GetVelGrad();
						dvxdx=velGrad->xx;
						dvydy=velGrad->yy;
						dvxdy=velGrad->xy;
						dvydx=velGrad->zz;			// yx stored in zz
						f2ForJx+=rho*((ax*duxdx+ay*duydx)-(vx*dvxdx+vy*dvydx)); 
						f2ForJy+=rho*((ax*duxdy+ay*duydy)-(vx*dvxdy+vy*dvydy));
					}
				}
				
				if(count==0)
					throw "J Integral contour contains no particles";
				carea=1.e-9*(cxmax-cxmin)*(cymax-cymin)/count;	// area per particle
				Jx2=f2ForJx*carea;      // Jx2 in Nmm/mm^2 now
				Jy2=f2ForJy*carea;      // Jy2 in Nmm/mm^2 now
			}
			
			/* Task 6: Subtract energy due to tractions or for cracks in contact, subtract
				energy associated with shear stress (not yet implemented thought)
				If Jterms==3 or 4, subtract recoverable energy as well for R curve analysis
			*/
			if(hasTractionLaws)
			{	tractionEnergy=startSeg->TractionEnergy(&crossPt,crkTipIdx,true);
				bridgingReleased=startSeg->TractionEnergy(&crossPt,crkTipIdx,false);
			}
			else
			{	// set traction energy to energy due to shear if in contact
				tractionEnergy=0.;
				bridgingReleased=0.;
			}
			
			// add the two terms N mm/mm^2
			Jx=Jx1+Jx2;
			Jy=Jy1+Jy2;

			/* Jint -- crack-axis components of dynamic J-integral
				  Jint.x is J1 in archiving and literature and is energy release rate, here
						it accounts for traction laws. Friction can be handled but not yet implemented
				  Jint.y literature J2 - needed only to convert to KI and KII (archive when propagation is off)
				  Jint.z is actual energy released when the crack and traction zone propagate together
							(archived as J2 when propagation is on)
			   crackDir -- crack propagating direction cosines from above
			*/
			tipCrk->Jint.x= Jx*crackDir.x+Jy*crackDir.y-tractionEnergy;		// Jtip or energy that will be released if crack grows
			tipCrk->Jint.y=-Jx*crackDir.y+Jy*crackDir.x;			// J2(x)
			//tipCrk->Jint.y= Jx*crackDir.x+Jy*crackDir.y;			// store J1(x) in traction zones (archived only when no propagation)
			tipCrk->Jint.z= tipCrk->Jint.x+bridgingReleased;		// Jrel or energy released in current state
			
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
	long iel;
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
		theElements[iel]->GetShapeFunctions(&numnds,fn,nds,&cpos,&cncpos);
	
        // Add crack particle heating to each node in the element
        for(i=1;i<=numnds;i++)
		{	nd[nds[i]]->fcond+=scrk->HeatRate()*fn[i];
		}
		
		// next segment
        scrk=scrk->nextSeg;
	}
}

// crack crossing method

// Determine if line from particle (1) to node (2) crosses this crack
// Return ABOVE_CRACK (1), BELOW_CRACK (2), or NO_CRACK (0) and crack normal in norm
short CrackHeader::CrackCross(double x1,double y1,double x2,double y2,Vector *norm)
{
    double x3,y3,x4,y4;
    CrackSegment *scrk=firstSeg;
    short cross=NO_CRACK;
    double area1,area2;
    
    // in no segments
    if(scrk==NULL) return cross;
    
    // check extents
    if(fmax(x1,x2)<xmin) return cross;
    if(fmin(x1,x2)>xmax) return cross;
    if(fmax(y1,y2)<ymin) return cross;
    if(fmin(y1,y2)>ymax) return cross;
    
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
            
        // first two areas (123 and 124)
        area1=Triangle(x1,y1,x2,y2,x3,y3);
        area2=Triangle(x1,y1,x2,y2,x4,y4);
            
        // check for crossing
        while(TRUE)
        {   // first area negative
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
}

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
