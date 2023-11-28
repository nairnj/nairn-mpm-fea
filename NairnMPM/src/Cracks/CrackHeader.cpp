/********************************************************************************
    CrackHeader.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Apr 3 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include <fstream>

#include "NairnMPM_Class/NairnMPM.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Cracks/CrackHeader.hpp"
#include "Materials/MaterialBase.hpp"
#include "System/ArchiveData.hpp"
#include "Elements/ElementBase.hpp"
#include "Cracks/CrackSurfaceContact.hpp"
#include "Cracks/ContourPoint.hpp"
#include "Cracks/CrackSegment.hpp"
#include "Cracks/CrackNode.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Cracks/CrackLeaf.hpp"
#include "Cracks/CrossedCrack.hpp"
#include "Read_XML/ParseController.hpp"
#include "Custom_Tasks/PropagateTask.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Exceptions/MPMWarnings.hpp"

typedef struct {
    int crackNum;
    int crackSign;
    CrackSegment *startSeg;
    NodalPoint *phantom;
} OtherCrack;


using namespace std; 

// class statics
double CrackHeader::codLocation;
int CrackHeader::codInterval;
double CrackHeader::bezArg[4];
double CrackHeader::bezDer[4];
int CrackHeader::warnNodeOnCrack;
int CrackHeader::warnThreeCracks;
int CrackHeader::warn2ndTipInContour;
int CrackHeader::warn2ndTipInCell;
int CrackHeader::warnAdatedContourFailed;

// globals
CrackHeader *firstCrack;		// first crack
int JGridSize = 2;				// size of J Integral contour
int JContourType = AXISYM_BROBERG_J;			// different methods in axisymmetric J integral
int JTerms = -1;				// number of terms in J Integral calculation (default 1 or 2 if axisymmetric)
int JGridEnergy = 0;			// Calculate work and kinetic energy on the grid GRID_JTERMS
int numberOfCracks = 0;
CrackHeader **crackList;

#pragma mark CrackHeader: Constructors and Destructors

// Constructors
CrackHeader::CrackHeader()
{	
    firstSeg = NULL;
    lastSeg = NULL;
    numberSegments = 0;
	numberKeypoints = -1;
	fixedCrack = false;
	customContactLawID = -1;
	customTractionPropID = -1;
	hasTractionLaws = false;
	thickness = 1.0;				// for crack tip heating and tractions in mm, will default to grid thickness if set
	allowAlternate[0] = allowAlternate[1] = true;
	crossedCracks = NULL;
	numberFacets = 0;				// only for 3D and if 0 it is 2D crack
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
// 1. Verify has segments
// 2. Check all segments for acceptable traction law (if used)
// 2a. If has traction laws, initialze law history data
// 2b. If any segments too shart, print warning
// 3. Check crack tip materials to be acceptable solid materials
// d. If shiftOption true, shift cracks slightly.
// throws CommonException()
// Note that 3D cracks calls this, so must segregate 2D only calcs.
void CrackHeader::PreliminaryCrackCalcs(double dcell,bool shiftOption,bool is3D)
{
    if(!is3D)
    {   // it does not make sense unless there are two segments (and at least one line)
        if(firstSeg==NULL)
        {	throw CommonException("A defined crack does not have any particles","CrackHeader::PreliminaryCrackCalcs");
            return;
        }
        else if(firstSeg->nextSeg==NULL)
        {	throw CommonException("All cracks must have at least two particles","CrackHeader::PreliminaryCrackCalcs");
            return;
        }
    }
    
	// loop through segments
	CrackSegment *scrk = firstSeg;
    double lastx = scrk->cp.x;
    double lasty = scrk->cp.y;
	bool tooShort = false, tooLong = false;
	double shortest,longest;
    
	while(scrk!=NULL)
	{   // check segment length
		if(scrk!=firstSeg && !is3D)
		{	double x = scrk->cp.x;
			double y = scrk->cp.y;
			double segLength = sqrt((x-lastx)*(x-lastx)+(y-lasty)*(y-lasty));
			if(segLength<0.1*dcell)
			{	if(tooShort)
					shortest=fmin(segLength,shortest);
				else
					shortest = segLength;
				tooShort = true;
			}
			if(segLength>2.*dcell)
			{	if(tooLong)
					longest=fmax(segLength,longest);
				else
					longest = segLength;
				tooLong = true;
			}
			lastx = x;
			lasty = y;
		}
		
		// check check traction laws and set history variables
		if(hasTractionLaws)
		{	int matid=scrk->MatID();
			if(matid>=0)
			{	if(matid>=nmat)
					throw CommonException("Crack segment with an undefined traction law material","CrackHeader::PreliminaryCrackCalcs");
				if(fixedCrack)
					throw CommonException("Fixed crack cannot have a traction law segment","CrackHeader::PreliminaryCrackCalcs");
				if(theMaterials[matid]->MaterialStyle()!=TRACTION_MAT)
					throw CommonException("Crack segment with material that is not a traction law","CrackHeader::PreliminaryCrackCalcs");
				
				// allow traction law to have history dependent data
				scrk->SetHistoryData(theMaterials[matid]->InitHistoryData(NULL,NULL));
			}
		}
		
		// next segment
		scrk=scrk->nextSeg;
	}
	
	// warning if too short or too long
	if(tooShort || tooLong)
	{	cout << "\nWARNING: At least one crack segment is ";
		if(tooShort) cout << "less than 1/10 a cell (shortest=" << shortest << ")";
		if(tooLong)
		{	if(tooShort) cout << " and one is ";
			cout << "longer than 2 cells (longest=" << longest << ")";
		}
		cout << "." << endl;
		cout << "         The simulation will run, but crack calculations may be inefficient" << endl;;
		cout << "         or inaccurate and crack CODs may be invalid." << endl;
	}
		
	// check crack tip materials to be valid and to not be a traction law material
    if(!is3D)
    {   CrackSegment *tipCrk;
        Vector crackDir;
        int crkTipIdx=START_OF_CRACK;
        while(crkTipIdx<=END_OF_CRACK)
        {	CrackTipAndDirection(crkTipIdx,&tipCrk,crackDir);
            if(tipCrk->tipMatnum>0)
            {	if(tipCrk->tipMatnum>nmat)
                    throw CommonException("Crack tip material is an undefined material","CrackHeader::PreliminaryCrackCalcs");
                if(theMaterials[tipCrk->tipMatnum-1]->MaterialStyle()!=SOLID_MAT)
                    throw CommonException("Crack tip material must be a solid material","CrackHeader::PreliminaryCrackCalcs");
            }
            crkTipIdx++;
        }
    }
}

#pragma mark CrackHeader: Set up Cracks

// add new crack segment to end of the list when creating cracks
// but if same as previous location (e.g., when connecting shapes), no
// need to add this one
bool CrackHeader::addSegment(CrackSegment *cs,bool ignoreDuplicate)
{
    if(cs==NULL) return false;		// not created
    if(lastSeg==NULL)
    {	firstSeg=cs;
    }
    else
	{	// no need to add a zero length segment (only checked in 2D)
		if(!ignoreDuplicate && !fmobj->IsThreeD())
		{	if(DbleEqual(cs->cp.x,lastSeg->cp.x) && DbleEqual(cs->cp.y,lastSeg->cp.y))
			{	// but maybe want new tip material or traction law
				if(cs->tipMatnum>0)
					lastSeg->tipMatnum=cs->tipMatnum;
				if(cs->MatID()>=0)
				{   lastSeg->SetMatID(cs->MatID()+1);
					hasTractionLaws = true;
				}
				return true;
			}
		}
    	lastSeg->nextSeg=cs;
		cs->prevSeg=lastSeg;
    }
    lastSeg=cs;
    numberSegments++;
	
    // determine element for new crack particles
    cs->FindInitialElement();
	if(cs->planeElemID()<0) return false;
	
	// has it put traction laws on this crack
	if(cs->MatID()>=0) hasTractionLaws = true;
    
    // give segment chance to store its hader
    cs->SetCrackHeader(this);
	
   return true;
}

// Change segment cohesive zone when meshing 3D cracks
void CrackHeader::setSegmentCZM(int czmID,int num,int numSet)
{
	// exit if not a segment number
	if(firstSeg==NULL || num>=numberSegments) return;
	
	// find segment number num
	CrackSegment *changeSeg = firstSeg;
	while(num>0)
	{	changeSeg = changeSeg->nextSeg;
		num--;
	}
	
	// change numSet segments starting with changeSeg (but not beyond lastSeg)
	while(numSet>0 && changeSeg!=NULL)
	{	changeSeg->SetMatID(czmID);
		changeSeg = changeSeg->nextSeg;
		numSet--;
	}
}

// output crack info and evaluate contact law
void CrackHeader::Output(void)
{
	cout << "  Crack " << number;
	if(IsThreeD())
		cout << ": points = " << NumberOfSegments() << ", facets = " << NumberOfFacets();
	else
		cout << ": length = " << Length() << ", segments = " << NumberOfSegments() << ", thickness = " << thickness << " mm";
	if(fixedCrack) cout << " (fixed)";
	if(hasTractionLaws)
	{	cout << " (has traction laws)";
		// if has traction laws, must use frictionless contact
		// When materials array created, last one (1 based) set to needed frictionless law if any traction laws

		customContactLawID = nmat;
	}
	cout << endl;
	
	// Crack tip settings in 2D
	if(!IsThreeD())
	{	if(firstSeg->tipMatnum!=-1)
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
		
		// save initial crack tip directions
		CrackSegment *crkTip;
		CrackTipAndDirection(START_OF_CRACK,&crkTip,initialDirection[START_OF_CRACK]);
		CrackTipAndDirection(END_OF_CRACK,&crkTip,initialDirection[END_OF_CRACK]);
		
		// check thickness
		double gridThickness=mpmgrid.GetThickness();
		if(gridThickness>0. && !DbleEqual(gridThickness,thickness))
			cout << "     +++++ WARNING: crack thickness does not match grid thickness +++++" << endl;
	}
	
	// custom contact law
	contact.CustomCrackContactOutput(customContactLawID,number);
	
	// custom traction propagation material
	if(customTractionPropID>0)
	{	// check OK
		if(customTractionPropID>nmat)
			throw CommonException("Propagation traction law for a crack is not defined","CrackHeader::Output");
		if(theMaterials[customTractionPropID-1]->MaterialStyle()!=TRACTION_MAT)
			throw CommonException("Propagation traction law for a crack is not a traction law material","CrackHeader::Output");
		cout << "    Custom propagation traction: " << theMaterials[customTractionPropID-1]->name
				<< " (number " << customTractionPropID << ")" << endl;
	}
}

#pragma mark CrackHeader: 2D-only methods

/* add new crack segment for crack propagation
	whichTip=0 (start) or 1 (for end)
    
    1. Add new segment at start or end
    2. Adjust firstSeg or lastSeg and nextSeg near the end points
    3. Set tipMatnum new segment to previous tip and set previous one to -1
    4. Transfer Crack tip properties to new crack tip
    5. Locate element for new crack tip
*/
bool CrackHeader::addSegmentTip(CrackSegment *cs,int whichTip)
{
    CrackSegment *prevSeg;
    double dx,dy;
    
    if(cs==NULL) return false;		// not created
	
	// find the element first (problem if not in the mesh)
    cs->FindInitialElement();
	
	// if it is not in the mesh, remove the segment and stop propagation
	if(cs->planeElemID()<0)
	{	if(whichTip==END_OF_CRACK)
			lastSeg->tipMatnum=-1;
        else
			firstSeg->tipMatnum=-1;
		delete cs;
		return false;
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
    
    // calculate growth that formed this segment
    dx=cs->cp.x-prevSeg->cp.x;
    dy=cs->cp.y-prevSeg->cp.y;
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
        cs->SetHistoryData(theMaterials[tmatnum]->InitHistoryData(NULL,NULL));
    }
	
    numberSegments++;
    
    ExtendHierarchy(cs);
    
    return true;
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
	{	cod.x=crkTip->surf[0].x-crkTip->surf[1].x;
		cod.y=crkTip->surf[0].y-crkTip->surf[1].y;
	}
	else
	{	cod.x=adjTip->surf[0].x-adjTip->surf[1].x;
		cod.y=adjTip->surf[0].y-adjTip->surf[1].y;
	}
	if(whichTip==START_OF_CRACK)
	{	cod.x=-cod.x;
		cod.y=-cod.y;
	}
	
	// convert cod.x to shear (mode II) cod and cod.y to normal (mode I)
	if(getModes)
	{	CrackSegment *newTip;
		CrackTipAndDirection(whichTip,&newTip,tipDir);
		double codx=cod.x, cody=cod.y;
		cod.x = codx*tipDir.x + cody*tipDir.y;
		cod.y = -codx*tipDir.y + cody*tipDir.x;
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
		{	cod.x=crkTip->surf[0].x-crkTip->surf[1].x;
			cod.y=crkTip->surf[0].y-crkTip->surf[1].y;
		}
		else
		{	cod.x=seg[1]->surf[0].x-seg[1]->surf[1].x;
			cod.y=seg[1]->surf[0].y-seg[1]->surf[1].y;
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
			cod.x=codx*tipDir.x + cody*tipDir.y;
			cod.y=-codx*tipDir.y + cody*tipDir.x;
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
		double codx=cod.x, cody=cod.y;
		cod.x = codx*tipDir.x + cody*tipDir.y;
		cod.y = -codx*tipDir.y + cody*tipDir.x;
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
		{	spt[i].x=seg[i]->cp.x;
			spt[i].y=seg[i]->cp.y;
		}
		else
		{	spt[i].x=seg[i]->surf[j].x;
			spt[i].y=seg[i]->surf[j].y;
		}
	}
	
	// calculate control points
	bpt[0].x=spt[0].x;
	bpt[0].y=spt[0].y;
	bpt[1].x=(-4.*spt[0].x+24.*spt[1].x-6.*spt[2].x+spt[3].x)/15.;
	bpt[1].y=(-4.*spt[0].y+24.*spt[1].y-6.*spt[2].y+spt[3].y)/15.;
	bpt[2].x=(spt[0].x-6.*spt[1].x+24.*spt[2].x-4.*spt[3].x)/15.;
	bpt[2].y=(spt[0].y-6.*spt[1].y+24.*spt[2].y-4.*spt[3].y)/15.;
	bpt[3].x=spt[3].x;
	bpt[3].y=spt[3].y;
	
#else
	
	// for bezier curves, crack particles are control points
	for(i=0;i<4;i++)
	{	if(surface==NO_CRACK)
		{	bpt[i].x=seg[i]->cp.x;
			bpt[i].y=seg[i]->cp.y;
		}
		else
		{	bpt[i].x=seg[i]->surf[j].x;
			bpt[i].y=seg[i]->surf[j].y;
		}
	}
	
	// calculate the interpolant points
	spt[0].x=bpt[0].x;
	spt[0].y=bpt[0].y;
	spt[1].x=(bpt[0].x+4.*bpt[1].x+bpt[2].x)/6.;
	spt[1].y=(bpt[0].y+4.*bpt[1].y+bpt[2].y)/6.;
	spt[2].x=(bpt[1].x+4.*bpt[2].x+bpt[3].x)/6.;
	spt[2].y=(bpt[1].y+4.*bpt[2].y+bpt[3].y)/6.;
	spt[3].x=bpt[3].x;
	spt[3].y=bpt[3].y;
	
#endif
	
	// control points for bezier curve between spt[codInterval] and spt[codInterval+1]
	// codInterval is 0 to 2 for which B spline to use with 0 being the one
	//   that starts at the crack tip (and the default setting)
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
	{	// default implementation has codInterval=0, bezArg = {0,0,0,1}
		// pos = bc[3] = spt[1] = (bpt[0]+4bpt[1]+bpt[2])/3
		// Default method sets bpt to crack points, so pos is weighted average of two ends of first two segment.
		for(i=0;i<4;i++)
		{	pos.x+=bezArg[i]*bc[i].x;
			pos.y+=bezArg[i]*bc[i].y;
		}
	}
}

// load vector with initial crack tip direction
void CrackHeader::GetInitialDirection(CrackSegment *crkTip,Vector &tipDir) const
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
	dx=(*crkTip)->cp.x-adjCrk->cp.x;
	dy=(*crkTip)->cp.y-adjCrk->cp.y;
#else

#ifndef _CUBIC_INTERPOLATION_

	// Bezier curve is tangent to last segment, thus direction is same as last segment
	// This is method that is currently active
	if(crkTipIdx==START_OF_CRACK)
	{	*crkTip=firstSeg;
		adjCrk=firstSeg->nextSeg;
	}
	else
	{	*crkTip=lastSeg;
		adjCrk=lastSeg->prevSeg;
	}
	dx=(*crkTip)->cp.x-adjCrk->cp.x;
	dy=(*crkTip)->cp.y-adjCrk->cp.y;
	
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
		dx+=19.*(adjCrk->cp.x-nextCrk->cp.x);			// S0-S1
		dy+=19.*(adjCrk->cp.y-nextCrk->cp.y);
		
		adjCrk=nextCrk->nextSeg;				// S2
		if(adjCrk!=NULL)
		{	dx+=5.*(adjCrk->cp.x-nextCrk->cp.x);		// S2-S1
			dy+=5.*(adjCrk->cp.y-nextCrk->cp.y);
			
			nextCrk=adjCrk->nextSeg;			// S3
			if(nextCrk!=NULL)
			{	dx+=adjCrk->cp.x-nextCrk->cp.x;		// S2-S3
				dy+=adjCrk->cp.y-nextCrk->cp.y;
			}
		}
	}
	else
	{	*crkTip=lastSeg;
	
		adjCrk=lastSeg;						// S0
		nextCrk=lastSeg->prevSeg;				// S1
		dx+=19.*(adjCrk->cp.x-nextCrk->cp.x);			// S0-S1
		dy+=19.*(adjCrk->cp.y-nextCrk->cp.y);
		
		adjCrk=nextCrk->prevSeg;				// S2
		if(adjCrk!=NULL)
		{	dx+=5.*(adjCrk->cp.x-nextCrk->cp.x);		// S2-S1
			dy+=5.*(adjCrk->cp.y-nextCrk->cp.y);
			
			nextCrk=adjCrk->prevSeg;			// S3
			if(nextCrk!=NULL)
			{	dx+=adjCrk->cp.x-nextCrk->cp.x;		// S2-S3
				dy+=adjCrk->cp.y-nextCrk->cp.y;
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

// The crack should propage to (xnew,ynew)
//	whichTip= 0 or 1 for start or end of crack
// throws std::bad_alloc
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
	propSeg=new CrackSegment(&grow,tipMatID,tractionMat);
	if(!addSegmentTip(propSeg,whichTip))
	{	cout << "# Crack No. " << number << " left the grid; propagation was stopped" << endl;
		return NULL;
	}
	if(tractionMat>0)
	{	hasTractionLaws = true;
		fmobj->SetHasTractionCracks(true);
	}
	
	// return the segment
	return propSeg;
}

// Add crack tip heating for all points in this crack
void CrackHeader::CrackTipHeating(void)
{
	CrackSegment *scrk=firstSeg;
	int iel;
	int i;
#ifdef CONST_ARRAYS
	double fn[MAX_SHAPE_NODES];
	int nds[MAX_SHAPE_NODES];
#else
	double fn[maxShapeNodes];
	int nds[maxShapeNodes];
#endif
	
	// exit if no segments
	if(scrk==NULL) return;
	
	// loop over crack tips
	while(scrk!=NULL)
	{	// get element and shape function to extrapolate to the node
		iel=scrk->planeElemID();
		theElements[iel]->GetShapeFunctionsForCracks(fn,nds,scrk->cp);
		int numnds = nds[0];
	
		// normalize shape functions
		double fnorm = 0.;
		for(i=1;i<=numnds;i++)
		{	if(nd[nds[i]]->NodeHasNonrigidParticles())
				fnorm += fn[i];
		}
		
		// Add crack particle heating to each node in the element
		for(i=1;i<=numnds;i++)
		{	conduction->AddFluxCondition(nd[nds[i]],scrk->HeatRate()*fn[i]/fnorm,false);
		}
		
		// next segment
		scrk=scrk->nextSeg;
	}
}

#pragma mark CrackHeader: Methods

// archive crack to file
// throws CommonException()
void CrackHeader::Archive(ofstream &afile)
{
    int i=0;
    CrackSegment *mseg=firstSeg;
    
    // create space for this crack
	int recSize=archiver->GetRecordSize();
    long blen=recSize;
	char *aptr = new (std::nothrow) char[blen];
    if(aptr==NULL)
        throw CommonException("Memory error writing crack data.","CrackHeader::Archive");
    
    // move to start and then draw
    bool threeD = fmobj->IsThreeD();
    while(mseg!=NULL)
    {	mseg->FillArchive(aptr,i,threeD);
	
 		// write this crack segment (ofstream should buffer for us)
		afile.write(aptr,blen);
		if(afile.bad())
			archiver->FileError("File error writing crack data","(results file)","CrackHeader::Archive");
			
		mseg=mseg->nextSeg;
        i++;
    }
    
	delete [] aptr;
 }

// If contact.GetMoveOnlySurfaces() is true, crack plane was already moved by
//		the surfaces; thus only need to move to midpoint and check if element has changed
// If contact.GetMoveOnlySurfaces() is false, move all crack plane particles
//		using CM velocities (precalculated and stored in field[0])
// Also must recalculate extent of crack in cnear[i] and cfar[i]
short CrackHeader::MoveCrack(void)
{
	CrackSegment *scrk=firstSeg;
	// move only surfaces
	if(contact.GetMoveOnlySurfaces())
	{	while(scrk != NULL)
		{
			if(!fixedCrack)
			{	// move to midpoint between upper and lower surface
				scrk->MovePosition();
				
				// did element move
				if(!scrk->FindElement()) return false;
				
				// make sure surface are on correct side of the crack
				if(contact.GetPreventPlaneCrosses())
				{	if(!scrk->CheckSurfaces()) return false;
				}
			}

			// next segments
			scrk = scrk->nextSeg;
		}
	}
	
	// move crack plane particles by CM velocity
	else
	{	int iel;
		double shapeNorm;
		int j,nodeCounter;
		Vector vcm,acm,delv,dela;
#ifdef CONST_ARRAYS
		double fn[MAX_SHAPE_NODES];
		int nds[MAX_SHAPE_NODES];
#else
		double fn[maxShapeNodes];
		int nds[maxShapeNodes];
#endif
		
		// loop over crack points
		while(scrk != NULL)
		{
			if(!fixedCrack)
			{	// get element and shape functions to extrapolate to the particle
				iel=scrk->planeElemID();
				theElements[iel]->GetShapeFunctionsForCracks(fn,nds,scrk->cp);
				int numnds = nds[0];
				
				// initialize
				ZeroVector(&delv);
				ZeroVector(&dela);
				nodeCounter=0;
				shapeNorm=0.;
				
				// extrapolate to crack particle
				for(j=1;j<=numnds;j++)
				{	if(nd[nds[j]]->GetCMVelocityTask8(&vcm,&acm))
					{	AddScaledVector(&delv,&vcm,fn[j]);
						AddScaledVector(&dela,&acm,fn[j]);
						nodeCounter++;
						shapeNorm+=fn[j];
					}
				}
				
				// move it or collapse it
				if(nodeCounter>0)
				{	// pass unnormalized velocity and acceleration to segment
					scrk->MovePosition(&delv,&dela,timestep,shapeNorm);
					
					// did element move? But, if leaves grid, we assume a round off and try
					// to revert to moving at the midplane of the two surfaces
					if(!scrk->FindElement())
					{	scrk->MovePositionToMidpoint();
						if(!scrk->FindElement()) return false;
					}
					
					// check surfaces
					if(contact.GetPreventPlaneCrosses())
					{	if(!scrk->CheckSurfaces()) return false;
					}
				}
				
				else if(contact.GetPreventPlaneCrosses())
				{	// crack in free space, but check if surfaces have moved
					if(!scrk->CheckSurfaces()) return FALSE;
				}
			}

			// next segments
			scrk=scrk->nextSeg;
		}
	}
    
    // move crack hierarchy
    MoveHierarchy();
	
    return TRUE;
}

// Move one crack surface according to current velocities
// side==ABOVE_CRACK (1) or BELOW_CRACK (2)
short CrackHeader::MoveCrack(short side)
{
    CrackSegment *scrk=firstSeg;
    int iel;
    short js=side-1,nodeCounter,j;
	double fnorm;
	Vector svelnp1,surfAcc;
#ifdef CONST_ARRAYS
	double fn[MAX_SHAPE_NODES];
	int nds[MAX_SHAPE_NODES];
#else
	double fn[maxShapeNodes];
	int nds[maxShapeNodes];
#endif
    
    // loop over crack points
    while(scrk!=NULL)
	{
		if(!fixedCrack)
		{	// get element and shape functions to extrapolate to the particle
			iel = scrk->surfaceElemID(side);			// now zero based
			theElements[iel]->GetShapeFunctionsForCracks(fn,nds,scrk->surf[js]);
			int numnds = nds[0];

			// get S v+ or PIC velocity in svelnp1
			ZeroVector(&svelnp1);
			ZeroVector(&surfAcc);
			fnorm = 0;
			nodeCounter = 0;
			
			// extrapolate those with velocity to the particle
			for(j=1;j<=numnds;j++)
			{	if(nd[nds[j]]->IncrementDelvSideTask8(side,number,fn[j],&svelnp1,&surfAcc,&fnorm,scrk))
					nodeCounter++;
			}
			
			// svelnp1 is Sum(fi pi) and fnorm = Sum(fi mi), normalize to velocity
			if(nodeCounter>0)
			{	ScaleVector(&svelnp1,1./fnorm);
                // surfAcc is Sum(fi fToti), XPIC(1) needs to normalize to an acceleration
				ScaleVector(&surfAcc,1./fnorm);
			}
			else
				nodeCounter=0;
			
			// Update crack surface velocity and position
			if(scrk->MoveSurfacePosition(side,&svelnp1,&surfAcc,timestep,(nodeCounter>0)))
			{   // if bottom moved above, check if above moved elements
                if(!scrk->FindElement(ABOVE_CRACK)) return false;
			}
			
			// did surface move elements
			if(!scrk->FindElement(side))
			{	return false;
			}
		}
            
        // on to next segement
        scrk=scrk->nextSeg;
    }
    return true;
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
void CrackHeader::AddTractionForce(void)
{   if(!hasTractionLaws) return;
    
	CrackSegment *cs=firstSeg;
	while(cs!=NULL)
	{	cs->AddTractionForceSeg(this);
		cs=cs->nextSeg;
	}
}

#pragma mark CrackHeader: J Integral 2D

// J-integral calculation (YJG)
void CrackHeader::JIntegral(void)
{
    int gridElem,gridNode,nextNearest,edgeDir,centerNode;
    int i,j;
    int JCSize[4];
    ContourPoint *crackPt,*prevPt,*nextPt,*firstOtherPt,*prevOtherPt;
    int crkTipIdx;
    CrackSegment *tipCrk,*excludedTip;
    double Jx,Jy,Jx1,Jy1,Jx2,Jy2,JxAS2;
	Vector crackDir;
	bool secondTry;
	
    // Clear previously cross cracks list. It is only used to make sure a propagation
    //      does go through another crack. If growth crosses a cross cracks, it
    //      is stopped on that crack in J says to keep going
    // Note: this linked list is for cracks crossed by the contours. It has cracks
    //      from both tips in one list (i.e., not cleared). The fact that
    //      both tips are in one list is not an issue, because the list will be short and
    //      checking will be efficient even if a lot of cracks in the simulation
    if(crossedCracks!=NULL) crossedCracks->ClearObjects();
    CrossedCrack *lastCrossedCrack = NULL;
        
	//---------------------------------------------------------
    // Calculate J-integrals for each crack tip (if activated)
	//---------------------------------------------------------
	
	// it may try two contours at each crack tip. First try is at NearestNode().
	// if that path crosses the crack twice, it tries a contour from the next nearest node.
	secondTry = false;
	
    // crack contour radii in +x,+y,-x,-y (i.e., edge 0 to 3 for d node element)
    // Future could input these more general contours
    /* hack to XFEM example to use 8 cell contours
    if(mtime>0.59e-3 && JGridSize<8)
    {   cout << "# Switch to 8 cell J contour" << endl;
        JGridSize=8;
    }
    */
    for(i=0;i<4;i++) JCSize[i]=JGridSize;
    excludedTip = NULL;

	// Look at each crack tip
	crkTipIdx=START_OF_CRACK;
	while(crkTipIdx<=END_OF_CRACK)
    {
#pragma mark J_Task_1:Find crack tip
		//---------------------------------------------------------
		// TASK 1: find crack tip and the crack direction
		//---------------------------------------------------------
		CrackTipAndDirection(crkTipIdx,&tipCrk,crackDir);
        
        // done if not a crack tip, go to next tip
        if(tipCrk->tipMatnum<0)
		{	ZeroVector(&tipCrk->Jint);
			crkTipIdx++;
            continue;
        }
		double crackr = tipCrk->cp.x;			// for axisymmetric cracks
		
		// block to catch problems
		crackPt = NULL;					// first contour pt in the J contour
        firstOtherPt = NULL;            // first phantom point with another crack
        prevOtherPt = NULL;             // previous phantom point with another crack
		try
		{   // in case this tip object is rejected, store the last crossed crack
            lastCrossedCrack = crossedCracks!=NULL ? (CrossedCrack *)crossedCracks->lastObject : NULL;
            
#pragma mark J_Task_2:Find ccw nodal points
			//-------------------------------------------------------------------
			// Task 2: find ccw nodal points JCSize[edge] from crack tip nearest
			//           nodal point. Find orientation of each line segment
			//-------------------------------------------------------------------
			
			// get two nearest nodes to crack tip point
			gridElem=tipCrk->planeElemID();
			centerNode=theElements[gridElem]->NearestNode(tipCrk->cp.x,tipCrk->cp.y,&nextNearest);
			if(secondTry) centerNode=nextNearest;
            gridNode = centerNode;

			// set material type based on material near the tip.
			int oldnum = tipCrk->tipMatnum;
			tipCrk->FindCrackTipMaterial(oldnum);
			if(tipCrk->tipMatnum != oldnum)
			{	cout << "# crack " << GetNumber() << "'s tip left material " << oldnum
						<< " and entered material " << tipCrk->tipMatnum << endl;
			}
			
			// step to the edge of the J Integral contour
			gridNode=theElements[gridElem]->NextNode(gridNode,&edgeDir);
			for(i=1;i<JCSize[edgeDir];i++)
			{   gridElem=theElements[gridElem]->Neighbor(gridNode);
				if(gridElem<0)
					throw "J integral contour at a crack tip does not fit in the grid";
				gridNode=theElements[gridElem]->NextNode(gridNode);
			}
			
			// Walk around the countour (0 and 4 are half edges)
            //        2
            //    ---------
            //   |         |
            // 3 |    *    | 1  * is the crack tip
            //   |         |
            //    ----|----
            //    4     0
			// Note that | is starting point, not the crack crossing
			
			// set point at | to contour point
			prevPt = crackPt = new ContourPoint(nd[gridNode]);
			
            edgeDir = (edgeDir+1) % 4;
			int numPts=JCSize[edgeDir];
            Rect crect;
            crect.xmin=9e99; crect.xmax=-9e99; crect.ymin=9e99; crect.ymax=-9e99;
			for(j=0;j<5;j++)
			{	// set next node to a countour point
			    gridNode=theElements[gridElem]->NextNode(gridNode);
				nextPt=new ContourPoint(nd[gridNode]);
				
				// link to previous point
				prevPt->SetNextPoint(nextPt);
				prevPt=nextPt;
				
				// find numPts-1 more contour points
				for(i=1;i<numPts;i++)
				{	// get next element that shares face with face starting on grieNode
					gridElem=theElements[gridElem]->Neighbor(gridNode);
					if(gridElem<0)
						throw "J integral contour at a crack tip does not fit in the grid";
					
					// set next contour point
					gridNode=theElements[gridElem]->NextNode(gridNode);
					nextPt=new ContourPoint(nd[gridNode]);
					
					// link to previous point
					prevPt->SetNextPoint(nextPt);
					prevPt=nextPt;
				}
				
				// check corners for extent of contour rectangle
                crect.xmin=fmin(crect.xmin,nd[gridNode]->x);
                crect.xmax=fmax(crect.xmax,nd[gridNode]->x);
                crect.ymin=fmin(crect.ymin,nd[gridNode]->y);
                crect.ymax=fmax(crect.ymax,nd[gridNode]->y);

                // Get next edge +/-x sums the two x's, edge +/-y uses two y's
                // and adjust to half edge when about to do j=4 edge
                edgeDir = (edgeDir+1) % 4;
                if(j==3)
                {   // uses first half of the edge (minus 1)
                    numPts = JCSize[(edgeDir+2)%4] - 1;
                }
                else
                    numPts = JCSize[edgeDir] + JCSize[(edgeDir+2)%4];
			}
			
			// connect end to start
			prevPt->SetNextPoint(crackPt);

#pragma mark J_Task_3:Find this and other cracks crossing the contour
			//-------------------------------------------------------------------
			// Task 3: find crack crossing point and interacting crack crossings
			//-------------------------------------------------------------------
			
			// Contour cross this crack at crossPt through firstSeg. crossContourPt is
			// inserted into the contour and labeled as a phantomNode
			Vector crossPt;
			CrackSegment *startSeg = NULL;
			ContourPoint *crossContourPt = NULL;
            int crossCount = 0;

			// loop over nodes in the contour
			nextPt = crackPt;				// start tracing the contour
			while(true)
			{   // J integral path is along the countour line from node1 to node2
				NodalPoint *node1 = nextPt->node;
				NodalPoint *node2 = nextPt->nextPoint->node;
				
				// Check if line from node1 to node2 crosses this crack (in [0]) or another crack (in [1])?
				Vector crossPt1[2];
				int crackNum[2];
				double fract[2];
				CrackSegment *foundSeg[2];
				crackNum[0] = 0;
				crackNum[1] = 0;
                CrackHeader *otherCrack = NULL;
                
                // Does a segment of this crack cross line betwewen node 1 and node 2
				// If yes, return crack segment and (x,y) location in crossPt1
				foundSeg[0] = ContourCrossCrack(nextPt,&crossPt1[0]);
				if(foundSeg[0] != NULL)
				{	// it has crossed this crack
					crossCount++;
					if(crossCount==2)
					{   // error unless crossPt1=crossPt (previous one), which implies endpoints for two adjacent
						// contour intervals, but otherwise an error
						if(!(DbleEqual(crossPt.x,crossPt1[0].x) && DbleEqual(crossPt.y,crossPt1[0].y)))
						{	if(secondTry)
							{   throw "Two different crossings between J-path and a crack";
							}
							else
							{	// try again from a different node
								throw "";
							}
						}
					}
					else if(crossCount>2)
					{   // only gets here if found endpoints before and now cannot be another matching endpoint
						throw "Two different crossings between J-path and a crack";
					}
					else
					{   // save crossing point on target crack and segment
						// crossing between nextPt and nextPt->nextPoint
						crossContourPt = nextPt;
						crossPt.x = crossPt1[0].x;
						crossPt.y = crossPt1[0].y;
						// find crack particle closer to the crack tipstart
						startSeg = crkTipIdx==END_OF_CRACK ? foundSeg[0]->nextSeg : foundSeg[0] ;
						crackNum[0] = GetNumber();
                        // fract of distance from node 1 to 2 to the crossing point
						fract[0] = nextPt->Fraction(crossPt);
					}
				}
				
                // look for an interacting crack
                // This stops at first crack intersecting the segment and thus would miss more than
                //   one other crack intersecting the same segment (which is unlikely in resolved modeling)
                otherCrack = firstCrack;
                while(otherCrack!=NULL)
                {	// only check different cracks
                    if(otherCrack == this)
                    {   otherCrack = (CrackHeader *)otherCrack->GetNextObject();
                        continue;
                    }
                    
                    // does a segment of nextCrack cross line from node 1 to 2?
                    foundSeg[1] = otherCrack->ContourCrossCrack(nextPt,&crossPt1[1]);
                    if(foundSeg[1] != NULL)
                    {	// found crack crossing
                        crackNum[1] = otherCrack->GetNumber();
                        fract[1] = nextPt->Fraction(crossPt1[1]);
                        
                        // keep list of crossed cracks
                        CrossedCrack *nextCross = NULL;
                        if(crossedCracks==NULL)
                            crossedCracks = new ParseController();
                        else
                        {	// check if this crack is already in the list
                            nextCross = (CrossedCrack *)crossedCracks->firstObject;
                            while(nextCross!=NULL)
                            {	if(nextCross->crack==otherCrack) break;
                                nextCross = (CrossedCrack *)nextCross->GetNextObject();
                            }
                        }
                        if(nextCross==NULL)
                        {	crossedCracks->AddObject(new CrossedCrack(otherCrack));
                        }
                        
                        // exit when find first other crossed cracks
                        break;
                    }
                    otherCrack = (CrackHeader *)otherCrack->GetNextObject();
                }
				
				// if found two crossing cracks, reorder if needed
				if(crackNum[0]!=0 && crackNum[1]!=0)
				{
                    if(fract[0]>fract[1])
					{	double tmp = fract[0];
						fract[0] = fract[1];
						fract[1] = tmp;
						Vector tmpPt = crossPt1[0];
						crossPt1[0] = crossPt1[1];
						crossPt1[1] = tmpPt;
						int tmpCn = crackNum[0];
						crackNum[0] = crackNum[1];
						crackNum[1] = tmpCn;
						CrackSegment *tmpSeg = foundSeg[0];
						foundSeg[0] = foundSeg[1];
						foundSeg[1] = tmpSeg;
					}
				}
				
				// If this contour segment crosses this crack and or another crack then must break into
				// two or three contour segments othersize use field [0] on both nodes.
				// This applies to crossing this crack or
                // at most one other crack (other crack might cross in two places)
				for(int cn=0;cn<2;cn++)
				{	if(crackNum[cn]==0) continue;
					
					// Create a phantom node at the crossing point
					NodalPoint *phantom = NodalPoint::Create2DNode(crackNum[cn],crossPt1[cn].x,crossPt1[cn].y);
					phantom->PrepareForFields();
                    
                    // New way to get crackSign
					
					// J integral will use phatom[0] up to a phantom contour point and phantom[1] starting from it by:
					//    Interpolate [0] from node1 and [1] or [2] that crosses crackNum from node2 to phantom [0]
					//    Interpolate [0] from node2 and [1] or [2] that crosses crackNum from node1 to phantom [1]
                    // crackSign=-1 means node2 above crack, node1 below crack, crack direction exiting
                    // crackSign=1 means node1 above crack and node 2 below, crack direction entering
                    int crackSign;
					phantom->Interpolate(node1,node2,fract[cn],crackNum[cn],&crackSign);
					
					// insert new ContourPoint with the phantom node between nextPt and nextPt->nextPoint
                    // node number of the phantom node is the crack number
					ContourPoint *insertPt=new ContourPoint(phantom);
					insertPt->SetPhantomNode(true);
					insertPt->SetNextPoint(nextPt->nextPoint);
					nextPt->SetNextPoint(insertPt);
                    
                    // if for another crack, get other crack properties for integrating
                    if(crackNum[cn]!=GetNumber())
                    {
                        // new way to get crack sign
                        if(PtInRect(&foundSeg[cn]->cp,&crect))
                        {   // implies ->nextSeg is out of contour
                            foundSeg[cn]=foundSeg[cn]->nextSeg;
                            crackSign = -1;
                            if(PtInRect(&foundSeg[cn]->cp,&crect))
                               cout << "# found unexpected segment inside contour" << endl;
                            /*
                            if(foundSeg[cn]->prevSeg==NULL)
                            {   // foundSeg must be crack tip at start of crack in contour
                                crackSign=-1;
                                foundSeg[cn]=foundSeg[cn]->nextSeg;
                            }
                            else if(!PtInRect(&foundSeg[cn]->prevSeg->cp,&crect))
                            {   crackSign=1;
                                foundSeg[cn]=foundSeg[cn]->prevSeg;
                            }
                            else
                            {   crackSign=-1;
                                foundSeg[cn]=foundSeg[cn]->nextSeg;
                            }
                            */
                        }
                        else
                        {   // implies ->nextSeg is in the contour
                            crackSign = 1;
                            if(!PtInRect(&foundSeg[cn]->nextSeg->cp,&crect))
                               cout << "# found unexpected segment outside contour" << endl;
                            /*
                            if(foundSeg[cn]->prevSeg==NULL)
                            {   // foundSeg is crack tip just outside the countour
                                crackSign=1;
                            }
                            else if(PtInRect(&foundSeg[cn]->prevSeg->cp,&crect))
                            {   crackSign=-1;
                            }
                            else
                            {   crackSign=1;
                            }
                            */
                        }
                        
                        // look for next segment outside, but kink will be inside
                        bool exitPtIsKink = false;
                        if(crackSign==-1)
                        {   if(foundSeg[cn]->nextSeg!=NULL)
                                exitPtIsKink = PtInRect(&foundSeg[cn]->nextSeg->cp,&crect);
                        }
                        else
                        {   if(foundSeg[cn]->prevSeg!=NULL)
                                exitPtIsKink = PtInRect(&foundSeg[cn]->prevSeg->cp,&crect);
                        }
                        
                        // if first other crack set first one now
                        if(firstOtherPt==NULL) firstOtherPt = insertPt;
                        
                        // set crack details, set previous to pount to this, and find if exits contour
                        insertPt->SetOtherCrack(otherCrack,foundSeg[cn],crackSign,prevOtherPt,&crect);
                        prevOtherPt = insertPt;
                            
                        // did this exit at a previous countour point
                        // search up to this point
                        if(insertPt->endSeg!=NULL)
                        {   ContourPoint *matchPt = firstOtherPt;

                            // step through prior crossing points to see if matching
                            while(matchPt!=insertPt)
                            {   if(matchPt->startSeg==insertPt->endSeg && matchPt->crackSign!=insertPt->crackSign)
                                {   // this new point matches previous one
                                    // record it and turn off integration at this point
                                    matchPt->endingPt = insertPt;
                                    insertPt->crackSign = 0;
                                    break;
                                }
                                matchPt = matchPt->nextOther;
                            }
                        }
                    }
                    
                    // reset nextPt to new insertPt
                    nextPt = insertPt;
				}
				
				// on to next segment
				nextPt=nextPt->nextPoint;
				if(nextPt==crackPt) break;
			}
			
            // if startSeg is still NULL, then did not find any crossings
            if(startSeg == NULL)
            {    throw "A crack does not cross its J path";
            }
 
#pragma mark J_Task_4:Integrate along contour
			//-------------------------------------------------------------------
 			// Task 4: Loop over all segments and evaluate J integral (from the primary term)
			// Transform to crack plane and save results
			//-------------------------------------------------------------------
			DispField *sfld1,*sfld2;

			// Initialize contour integration terms
			Jx1 = Jy1 = 0.0;				// J-integral components from the first term
			double tractionEnergy = 0.,bridgingReleased = 0.;		// for traction laws
			nextPt = crackPt;				// start tracing the contour
			int count = 0;					// particles in the nodal fields
			double r1 = 1.,r2 = 1.;			// for axisymmetric J's

			// loop over nodes in the contour
			while(true)
			{   // J integral along the countour line from node1 to node2
				NodalPoint *node1 = nextPt->node;
				NodalPoint *node2 = nextPt->nextPoint->node;
				
				// Get fields
                int cloc;
				count += node1->GetFieldForCrack(nextPt->phantomNode,true,&sfld1,0,cloc);
				count += node2->GetFieldForCrack(nextPt->nextPoint->phantomNode,false,&sfld2,0,cloc);

				// get r for Bergkvist and Huong axisymmetric calcs (JContourType always AXISYM_BROBERG_J when not axisymmetric)
				if(JContourType != AXISYM_BROBERG_J)
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
				double fForJx1,fForJy1,fForJx2,fForJy2,Jxs,Jys;

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

				// term (ti*ui,x) (uN/mm^2)
				termForJx1=(sxx1*segNorm.x+sxy1*segNorm.y)*dudx1
						  +(sxy1*segNorm.x+syy1*segNorm.y)*dvdx1;
				termForJx2=(sxx2*segNorm.x+sxy2*segNorm.y)*dudx2
						  +(sxy2*segNorm.x+syy2*segNorm.y)*dvdx2;
						  
				// [(W+K)nx-ti*ui,x] (uN/mm^2)
				fForJx1=(wd1+kd1)*segNorm.x-termForJx1;
				fForJx2=(wd2+kd2)*segNorm.x-termForJx2;

				// add for two endpoints using midpoint rule (r1=r2=1 unless axisymmetry by Bergkvist)
				Jxs = 0.5*(r1*fForJx1 + r2*fForJx2)*ds;	// uN mm/mm^2
				Jx1 += Jxs;

				// calculate Jy (or Jz if axisymmetric)

				// term ti*ui,y (uN/mm^2)
				termForJy1 = (sxx1*segNorm.x+sxy1*segNorm.y)*dudy1
                            +(sxy1*segNorm.x+syy1*segNorm.y)*dvdy1;
				termForJy2 = (sxx2*segNorm.x+sxy2*segNorm.y)*dudy2
                            +(sxy2*segNorm.x+syy2*segNorm.y)*dvdy2;
						  
				// [(W+K)ny-ti*ui,y] (uN/mm^2)
				fForJy1=(wd1+kd1)*segNorm.y-termForJy1;
				fForJy2=(wd2+kd2)*segNorm.y-termForJy2;

				// add for two endpoints using midpoint rule (r1=r2=1 unless axisymmetry by Bergkvist)
				Jys = 0.5*(r1*fForJy1 + r2*fForJy2)*ds;	// uN mm/mm^2
				Jy1 += Jys;

				// on to next segment
				nextPt=nextPt->nextPoint;
				if(nextPt==crackPt) break;
			}
			
			// if did not find any velocity fields in the countour, then the entire countour is
			// in empty space and J will be zero (and meaningless)
			if(count==0)
			{	throw "Section of the J Integral contour was in empty space";
			}

#pragma mark J_Task_4b:Integrate along intersecting cracks
            // Code Hack: If following line is uncommented, calculations for
            // intersecting cracks will not be done. Comment out to include effects
            // of intersecting crack (only uncomment for testing). It also turns
            // off traction correction below
            //firstOtherPt = NULL;

            // 1. firstOtherPt==NULL, no crossings so skip this section
            // 2. firstOtherPt!=NULL, integrate crossing cracks until tip
            //      or until leaves the contour at another crossing
            // 3. If find on that goes to a tip try to adapt the contour
            //      to omit that tip (and then start over)
            prevOtherPt = firstOtherPt;
            while(prevOtherPt!=NULL)
            {   // skip those not being integrated
                if(prevOtherPt->crackSign==0)
                {   prevOtherPt = prevOtherPt->nextOther;
                    continue;
                }
                
                // Get J terms up and back along this crack
                Vector Jcrack = JIntersectionCrack(prevOtherPt);
                Jx1 += Jcrack.x;
                Jy1 += Jcrack.y;
                
                // If was a tip segment, try to adapt the contour intead
                // If possible, then start over, otherwise warning and proceed
                if(prevOtherPt->endSeg==NULL)
                {   // find crack tip segment for this intersecting crack
                    CrackSegment *scrk2 = prevOtherPt->startSeg;
                    CrackSegment *scrk1 = NULL;
                    while(scrk2!=NULL)
                    {   scrk1 = scrk2;
                        scrk2 = prevOtherPt->crackSign>0 ? scrk2->nextSeg : scrk2->prevSeg;
                    }
                    
                    // scrk1 is other crack tip; find spacing between tips
                    int otherTipElem = scrk1->planeElemID();
                    int orow,ocol,rank,trow,tcol;
                    mpmgrid.GetElementsRCR(otherTipElem,orow,ocol,rank);
                    mpmgrid.GetElementsRCR(tipCrk->planeElemID(),trow,tcol,rank);
                    int dx = ocol-tcol;
                    int dy = orow-trow;
                    
                    // issue warning first time two tips are found
                    if(warnings.Issue(CrackHeader::warn2ndTipInContour,0)==GAVE_WARNING)
                    {   Describe();
                        cout << "#    Element spacings between two tips dx=" << dx << " dy=" << dy << endl;
                    }
                    
                    if(dx==0 && dy==0)
                    {   // two tips in the same element
                        // warn first time but go ahead with current J
                        if(warnings.Issue(CrackHeader::warn2ndTipInCell,0)==GAVE_WARNING)
                        {   Describe();
                        }
                    }
                    else
                    {   // reduce the contour size in one direction
                        // check target node vs crack to make sure avoids other crack
                        // if possible increase distance from 2nd crack but do not
                        //      get close to first than 2 cells
                        if(abs(dx)>=abs(dy))
                        {   if(dx>0)
                            {   if(nd[centerNode]->x>tipCrk->cp.x) dx--;
                                if(dx>2) dx--;
                                if(dx>2) dx--;
                                JCSize[0]=dx;
                            }
                            else
                            {   if(nd[centerNode]->x<tipCrk->cp.x) dx++;
                                if(dx<-2) dx++;
                                if(dx<-2) dx++;
                                JCSize[2]=abs(dx);
                            }
                        }
                        else
                        {   // reduce vertical direction
                            if(dy>0)
                            {   if(nd[centerNode]->y>tipCrk->cp.y) dy--;
                                if(dy>2) dy--;
                                if(dy>2) dy--;
                                JCSize[1]=dy;
                            }
                            else
                            {   if(nd[centerNode]->y<tipCrk->cp.y) dy++;
                                if(dy<-2) dy++;
                                if(dy<-2) dy++;
                                JCSize[3]=abs(dy);
                            }
                        }
                        
                        if(scrk1==excludedTip)
                        {   // same tip still inside the zone
                            // should not happen if above algorithm is correct
                            // but could happen if both tips of one crack in the contour
                            // give warning but proceed with calculated J terms
                            if(warnings.Issue(CrackHeader::warnAdatedContourFailed,0)==GAVE_WARNING)
                            {   Describe();
                            }
                        }
                        else
                        {   // comment out to test non-adative contours
                            excludedTip = scrk1;
                            throw "C";
                        }
                    }
                }   // end block when intersecting crack has tip in contour
                    
                // on to next intersecting crack
                prevOtherPt = prevOtherPt->nextOther;
            }   // end block intergrating other cracks
            
            // If gets here, current contour has not other crack tips
            // or adapting the countour failue and we give up ain continue
            // with that contour
            excludedTip = NULL;

#pragma mark J_Task_5:Evaluate J integral from the additional terms (dynamic and axisymmetric)
			Jx2 = Jy2 = JxAS2 = 0.;
			if(JTerms==2)
			{   double xp,yp,carea;
				double ax,ay,duxdx,duydx,duxdy,duydy;
				double vx,vy,dvxdx,dvydy,dvxdy,dvydx;
				double f2ForJx=0.,f2ForJy=0.,f2axisym=0.;
				double rho0,Jp,mp;
				double Vtot = 0.;	// volume of particles in the zone
				double Atot = 0.;	// area of particle for axisummetric term

				// integrate nonrigid particles
				for(int p=0;p<nmpmsNR;p++)
				{	MPMBase *mptr = mpm[p];
					if(mptr->InReservoir()) continue;
					xp=mptr->pos.x;
					yp=mptr->pos.y;
					if(XYInRect(xp,yp,&crect))
					{	// material reference
						const MaterialBase *matref = theMaterials[mptr->MatID()];
                        
						// if this material ignores cracks, then do not include it in J integral
						if(!matref->AllowsCracks()) continue;
						
						// get volume data and add to volume
						rho0 = matref->GetRho(mptr);
						Jp = matref->GetCurrentRelativeVolume(mptr,0);
						mp = mptr->mp;
						Vtot += Jp*mp/rho0;			// add particle volume
						
						// Accelerations mm/sec^2
						Vector *acc=mptr->GetAcc();
						ax=acc->x;
						ay=acc->y;
						
						// Displacement gradients (dimensionless)
						//Matrix3 gradU = mptr->GetDisplacementGradientMatrix();
						Matrix3 gradU = mptr->GetDisplacementGradientForJ(matref);
						duxdx = gradU(0,0);
						duydy = gradU(1,1);
						duxdy = gradU(0,1);
						duydx = gradU(1,0);
						
						// Velocities (mm/sec)
						vx=mptr->vel.x;
						vy=mptr->vel.y;
						
						// Velocity gradients (1/sec)
						Tensor *velGrad=mptr->GetVelGrad();
						dvxdx=velGrad->xx;
						dvydy=velGrad->yy;
						dvxdy=velGrad->xy;
						dvydx=velGrad->zz;			// yx stored in zz
						
						// Add an area integral is axisymetric
						if(fmobj->IsAxisymmetric())
						{	// Theory has area integral so need new volume element and track separate area
							double Ap = mptr->GetVolume(DEFORMED_AREA_FOR_GRADIENT);
							Atot += Ap;
							
							// in axisymmetrix z is theta direction, etheta = u/r. but w=0, az=vz=0
							// Since w=0, no change to above terms, but have some static terms for Jx=Jr only
							// Stress Units uN/mm^2 (mm^3/g)
							Tensor sp = mptr->ReadStressTensor();
							
							// Units uN/mm = uN mm/mm^2
							// add rho*() weighted by area or Ap*(rho0/Jp)*()
							if(JContourType == AXISYM_BROBERG_J)
							{	// See Broberg, Cracks and Fracture (1999), page 65
								f2axisym += (Ap*rho0/Jp)*(sp.xx*duxdx - sp.zz*gradU(2,2) + sp.xy*duydx)/xp;
								
								// We want Vp*rho*() = mp*()
								// increment the integrands g*(mm/sec^2) = uN
								f2ForJx += mp*((ax*duxdx+ay*duydx)-(vx*dvxdx+vy*dvydx));
								f2ForJy += mp*((ax*duxdy+ay*duydy)-(vx*dvxdy+vy*dvydy));
							}
							else
							{	// Bergkvist and Huong called J3D/(a dphi)
								f2axisym += (Ap*rho0/Jp)*(mptr->GetWorkEnergy() - sp.zz*gradU(2,2))/crackr;
								
								// We want Vp*rho*() = mp*()
								// increment the integrands g*(mm/sec^2) = uN
								f2ForJx += mptr->pos.x*mp*((ax*duxdx+ay*duydx)-(vx*dvxdx+vy*dvydx))/crackr;
								f2ForJy += mptr->pos.x*mp*((ax*duxdy+ay*duydy)-(vx*dvxdy+vy*dvydy))/crackr;
							}
						}
						else
						{	// We want Vp*rho*() = mp*()
							// increment the integrands g*(mm/sec^2) = uN
							f2ForJx += mp*((ax*duxdx+ay*duydx)-(vx*dvxdx+vy*dvydx));
							f2ForJy += mp*((ax*duxdy+ay*duydy)-(vx*dvxdy+vy*dvydy));
						}
					}
				}
				
				if(Vtot<=0.)
					throw "J Integral contour contains no particles";
                carea = RectArea(&crect);			            // Contour area mm^2
				Jx2 = f2ForJx*carea/Vtot;						// Jx2 in uN mm/mm^2 now
				Jy2 = f2ForJy*carea/Vtot;						// Jy2 in uN mm/mm^2 now
				if(fmobj->IsAxisymmetric())
					JxAS2 = f2axisym*carea/Atot;				// JxAS2 (for Jr in axisymmetric) in uN mm/mm^2
			}
			
#pragma mark J_Task_6:Crack tractions and contact
            /* Task 6: Subtract energy due to tractions or for cracks in contact, subtract
				energy associated with shear stress (later not yet implemented thought)
			*/
			if(hasTractionLaws)
			{	CrackSegment *tipSegment = startSeg;
				
				// Int(T.du) at current deformation
                // Warning: this calculation would be incorrect if segments with no traction law
                //   is between segments with traction law (withing the J contour)
				tractionEnergy = startSeg->TractionEnergy(&crossPt,crkTipIdx,true,&tipSegment);
                
                // look for intersections
                // Warning: this correction would be incorrect if there is a disruption of cohesive
                //  zone within the contour. It might also have small error if a crack crosses between
                //  tipSegment and edge of the J contour.
                //firstOtherPt = NULL;        // uncomment to compare to skipping these corrections
                if(firstOtherPt!=NULL)
                {   CrackSegment *prevSeg=tipSegment;
                    CrackSegment *closerSeg;
                    while(true)
                    {   closerSeg = (crkTipIdx==START_OF_CRACK) ? prevSeg->prevSeg : prevSeg->nextSeg ;
                        if(closerSeg==NULL) break;
                        
                        // check if line between segments is crossed by a crack (need to check all)
                        // if many cracks, could mark cracks as in the contour
                        CrackHeader *otherCrack = firstCrack;
                        Vector norm;
                        while(otherCrack!=NULL)
                        {    // only check different cracks
                            if(otherCrack == this)
                            {   otherCrack = (CrackHeader *)otherCrack->GetNextObject();
                                continue;
                            }
                            
                            // does this segment cross another crack
                            if(otherCrack->CrackCross(&prevSeg->cp,&closerSeg->cp,&norm,0)!=NO_CRACK)
                            {   tractionEnergy -= prevSeg->TractionEnergy(NULL,crkTipIdx,true,&tipSegment);
                                tractionEnergy += closerSeg->TractionEnergy(NULL,crkTipIdx,true,&tipSegment);
                                break;
                            }
                            
                            otherCrack = (CrackHeader *)otherCrack->GetNextObject();
                        }
                        
                        // next one
                        prevSeg = closerSeg;
                    }
                }
				
				// amount released so far in cohesive zone
				bridgingReleased = startSeg->TractionEnergy(&crossPt,crkTipIdx,false,NULL);
				
				// extra traction correction if Bergkvist and Huong axisymmetric J integral only
                // note that non-axisymmetri is always AXISYM_BROBERG_J
				if(JContourType != AXISYM_BROBERG_J)
				{	tractionEnergy *= crossPt.x/crackr;		// scale contour point energy
					
					// integrate tractions energy up to the crack tip
					CrackSegment *closerSeg=tipSegment;
					Vector previousPt;
					if(tipSegment==startSeg)
					{	// start at the contour point
						previousPt = crossPt;
					}
					else
					{	// zone is within the J contour, so start at root of the zone
						previousPt = tipSegment->cp;
						closerSeg = (crkTipIdx==START_OF_CRACK) ? closerSeg->prevSeg : closerSeg->nextSeg ;
					}
					double extra = 0.;
					while(closerSeg!=NULL)
					{	double dr = closerSeg->cp.x-previousPt.x;
						previousPt = closerSeg->cp;
						extra += closerSeg->TractionEnergy(&previousPt,crkTipIdx,true,NULL)*dr;
						closerSeg = (crkTipIdx==START_OF_CRACK) ? closerSeg->prevSeg : closerSeg->nextSeg ;
					}
					tractionEnergy += extra/crackr;
				}
			}
			else
			{	// set traction energy to energy due to shear if in contact (not implemented yet)
				// may also need to accout for contact compression (but that may be automatic in J integral
				//		if crack fields stresses are  continuous in normal direction on crack plane)
				tractionEnergy=0.;
				bridgingReleased=0.;
			}
			
#pragma mark J_Task_6:Finish up and store results
            
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
			   Units uN/mm, multiply by 1e-3 to get N/m = J/m^2
			*/
			tipCrk->Jint.x = Jx*crackDir.x + Jy*crackDir.y - tractionEnergy;		// Jtip or energy that will be released if crack grows
			tipCrk->Jint.y =-Jx*crackDir.y + Jy*crackDir.x;                         // J2(x) - for growth normal to crack plane
			tipCrk->Jint.z = tipCrk->Jint.x + bridgingReleased;						// Jrel or energy released in current state
			
			// end of try block on J calculation
			secondTry = false;
		}
		catch(const char *msg)
		{	// throwing "" signals to try again with the next nearest node
            // throwing "C" signals to try again with new contour
			if(strlen(msg)==0)
			{	secondTry = true;
                
                // contour is released below before trying again
                
                // release crossed cracks added on this failed attempt
                if(lastCrossedCrack!=NULL) crossedCracks->ClearObjects(lastCrossedCrack);
			}
            else if(msg[0]=='C')
            {   // contour is released below before trying again
                
                // release crossed cracks added on this failed attempt
                if(lastCrossedCrack!=NULL) crossedCracks->ClearObjects(lastCrossedCrack);
            }
			else
			{	cout << "# Crack No. " << number << ": "<< msg << endl;
				cout << "# J calculation and propagation at tip " << crkTipIdx << " will stop" << endl;
				tipCrk->tipMatnum=-1;
				ZeroVector(&tipCrk->Jint);
			}
		}
		catch(std::bad_alloc&)
		{	throw "Memory error in CrackHeader::JIntegral()";
		}
		catch( ... )
		{	throw "Unknown exception in CrackHeader::JIntegral()";
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
		
        // advance to next crack tip unless trying again for this tip
		if(!secondTry && excludedTip==NULL)
        {   crkTipIdx++;
            for(i=0;i<4;i++) JCSize[i]=JGridSize;
        }
        
    } // end loop over crack tips
}

// J integral on sides of a crack (do not call if crackSign in zero)
Vector CrackHeader::JIntersectionCrack(ContourPoint *cpt)
{
#ifdef CONST_ARRAYS
    double fn[MAX_SHAPE_NODES];
    int nds[MAX_SHAPE_NODES];
#else
    double fn[maxShapeNodes];
    int nds[maxShapeNodes];
#endif
	
    // factor to account for direction
    double scale = (double)cpt->crackSign;
    
    // Final Jx and Jy
    Vector Jcrack;
    ZeroVector(&Jcrack);
    
    // two displacement fields and tractions
    DispField *df2=new DispField;
    DispField *df1=new DispField;
    Vector T1,T2;
    bool useT1,useT2;
    
    // ABOVE_CRACK=1 and BELOW_CRACK=2
    for(int side=ABOVE_CRACK;side<=BELOW_CRACK;side++)
    {   short js = side-1;       // index to crack surface info

        // loop will exit when reach segment ending on cpt->endSeg or reach a crack tip
        CrackSegment *scrk2 = cpt->startSeg;
        CrackSegment *scrk1 = NULL;
        while(scrk2!=NULL)
        {   // -------- Task 1: extrapolate displacement field on grid to crack particle
            //                  depends on which side is being done
            // get element and shape functions to extrapolate to the particle
            int iel = scrk2->surfaceElemID(side);            // now zero based
            theElements[iel]->GetShapeFunctionsForCracks(fn,nds,scrk2->surf[js]);
            int numnds = nds[0];

            // clear needed df terms
            ZeroVector(&df2->du);
            ZeroVector(&df2->dv);
            df2->kinetic=0.;
            df2->work=0.;
            ZeroTensor(&df2->stress);
            
            // extrapolate those with displacement field to the particle
            double fnorm = 0.;
            int nodeCounter = 0;
            for(int j=1;j<=numnds;j++)
            {   if(nd[nds[j]]->IncrementDispField(side,number,fn[j],df2,scrk2))
                {   nodeCounter++;
                    fnorm += fn[j];
                }
            }
            
            // df is Sum(fi dfi) and fnorm = Sum(fi), normalize in case some empty nodes
            // for this crack surface.
            if(nodeCounter<numnds && fnorm>0.)
            {   double dscale = 1./fnorm;
                ScaleVector(&df2->du,dscale);
                ScaleVector(&df2->dv,dscale);
                df2->kinetic *= dscale;
                df2->work *= dscale;
                ScaleTensor(&df2->stress,dscale);
           }
            
           // -------- Task 2: check for crack traction (a traction law or zero if opened)
            useT2 = scrk2->SurfaceTraction(this,&T2);
            
            // -------- Task 3: first segment just stores it and goes on to next segment
            // first one just gets initial displacement field
            if(scrk1==NULL)
            {   *df1 = *df2;
                useT1 = useT2;
                T1 = T2;
                scrk1 = scrk2;
                scrk2 = cpt->crackSign>0 ? scrk2->nextSeg : scrk2->prevSeg;
                continue;
            }
            
            // -------- Task 4: Add to Jx and Jy for this crack segment and side

            // get length and the normal
            // normal accounts for side, scale accounts for counter clockwise contour direction
            double dx = scrk2->surf[js].x - scrk1->surf[js].x;
            double dy = scrk2->surf[js].y - scrk1->surf[js].y;
            double ds = sqrt(dx*dx+dy*dy);
            double nx,ny;
            if(side==ABOVE_CRACK)
            {   nx = dy/ds;
                ny = -dx/ds;
            }
            else
            {   nx = -dy/ds;
                ny = dx/ds;
            }
            
            // get traction if needed (only closed cracks)
            if(!useT1)
            {   T1.x = df1->stress.xx*nx + df1->stress.xy*ny;
                T1.y = df1->stress.xy*nx + df1->stress.yy*ny;
            }
            if(!useT2)
            {   T2.x = df2->stress.xx*nx + df2->stress.xy*ny;
                T2.y = df2->stress.xy*nx + df2->stress.yy*ny;
            }

            // calculate Jx (or Jr if axisymmetric by Broberg only)
#
            // [(W+K)nx-ti*ui,x]
            double fForJx1 = (df1->work+df1->kinetic)*nx - (T1.x*df1->du.x + T1.y*df1->dv.x);
            double fForJx2 = (df2->work+df2->kinetic)*nx - (T2.x*df2->du.x + T2.y*df2->dv.x);

            // add for two endpoints using midpoint rule (does not support Bergkvist for axisymmetric)
            // special case first and last segment
            double Jxs,phi;
            if(scrk1==cpt->startSeg)
            {   double dcpx = scrk2->cp.x - scrk1->cp.x;
                double dcpy = scrk2->cp.y - scrk1->cp.y;
                double dprex = cpt->node->x - scrk1->cp.x;
                double dprey = cpt->node->y - scrk1->cp.y;
                double dcl = sqrt(dcpx*dcpx+dcpy*dcpy);
                double dpre = sqrt(dprex*dprex+dprey*dprey);
                phi = dpre/dcl;
                Jxs = 0.5*(fForJx1*(1.-phi)*(1.-phi) + fForJx2*(1.-phi*phi))*ds;
            }
            else if(scrk2==cpt->endSeg)
            {   // if endingpPt==NULL the crossing never found a match
                // This is rare and might be kink (crosses back) or some other erro
                if(cpt->endingPt!=NULL)
                {   double dcpx = scrk2->cp.x - scrk1->cp.x;
                    double dcpy = scrk2->cp.y - scrk1->cp.y;
                    double dprex = cpt->endingPt->node->x - scrk1->cp.x;
                    double dprey = cpt->endingPt->node->y - scrk1->cp.y;
                    double dcl = sqrt(dcpx*dcpx+dcpy*dcpy);
                    double dpre = sqrt(dprex*dprex+dprey*dprey);
                    phi = dpre/dcl;
                    Jxs = 0.5*(fForJx1*(2.-phi)*phi + fForJx2*phi*phi)*ds;
                }
                else
                {   // just take whole segment
                    Jxs = 0.5*(fForJx1 + fForJx2)*ds;
                }
            }
            else
            {   // segment inside the contour
                Jxs = 0.5*(fForJx1 + fForJx2)*ds;
            }
            Jcrack.x += scale*Jxs;

            // calculate Jy (or Jz if axisymmetric by Broberg only)

            // [(W+K)ny-ti*ui,y]
            double fForJy1 = (df1->work+df1->kinetic)*ny - (T1.x*df1->du.y + T1.y*df1->dv.y);
            double fForJy2 = (df2->work+df2->kinetic)*ny - (T2.x*df2->du.y + T2.y*df2->dv.y);

            // add for two endpoints using midpoint rule (does not support Bergkvist for axisymmetric)
            double Jys;
            if(scrk1==cpt->startSeg)
                Jys = 0.5*(fForJy1*(1.-phi)*(1.-phi) + fForJy2*(1.-phi*phi))*ds;
            else if(scrk2==cpt->endSeg)
                Jys = 0.5*(fForJy1*(2.-phi)*phi + fForJy2*phi*phi)*ds;
            else
                Jys = 0.5*(fForJy1 + fForJy2)*ds;
            Jcrack.y += scale*Jys;
            
            // done this side if at endSeg
            if(scrk2==cpt->endSeg) break;
            
           // on to next segment
            *df1 = *df2;
            useT1 = useT2;
            T1 = T2;
            scrk1 = scrk2;
            scrk2 = cpt->crackSign>0 ? scrk2->nextSeg : scrk2->prevSeg;
        }
        
        // repeat for the other side
    }
    
    // clean up
    delete df1;
    delete df2;
    
    // return results
    return Jcrack;
}

// start in segment after startSeg and search for next crack plane position that is
// outside the rectangle. If get to tip first, return NULL
// csign=-1 means to step prevSeg and csign=1 means step nextSeg
CrackSegment *CrackHeader::FindRectExit(CrackSegment *startSeg,int csign,Rect *rect)
{
    CrackSegment *nextSeg = csign>0 ? startSeg->nextSeg : startSeg->prevSeg ;
    while(nextSeg!=NULL)
    {   if(!PtInRect(&nextSeg->cp,rect)) return nextSeg;
        nextSeg = csign>0 ? nextSeg->nextSeg : nextSeg->prevSeg ;
    }
    return NULL;
}

#pragma mark CrackHeader: 2D Crack crossing and hierachy

// Determine if line from particle xp1 to node xp2 crosses this crack
// Return ABOVE_CRACK (1), BELOW_CRACK (2), or NO_CRACK (0) and crack normal in norm
// This method uses hierarchical crack in a binary tree for 2D cracks
short CrackHeader::CrackCross(Vector *xp1,Vector *xp2,Vector *norm,int theNode) const
{
	// recursive method to traverse tree hierarchy
	return CrackCrossLeaf(rootLeaf,xp1,xp2,norm,NO_CRACK);
}

// Recursive Method to process each leaf in hierarchical traversal
short CrackHeader::CrackCrossLeaf(CrackLeaf *leaf,Vector *xp1,Vector *xp2,Vector *norm,short cross) const
{
    // check extents this leaf, return current cross if not in this leaf's extent
	if(!LineIsInExtents(xp1,xp2,leaf->cnear,leaf->cfar)) return cross;
    
    // It is in extent of this leaf
    // if not terminal, go on to the child leaves
    if(!leaf->ChildrenAreSegments())
    {   CrackLeaf *child1,*child2;
        leaf->GetChildLeaves(&child1,&child2);
        cross = CrackCrossLeaf(child1,xp1,xp2,norm,cross);
        if(child2!=NULL) cross = CrackCrossLeaf(child2,xp1,xp2,norm,cross);
        return cross;
    }
    
    // Terminal Leaf - check each segment
    CrackSegment *scrk1,*scrk2;
    leaf->GetChildSegments(&scrk1,&scrk2);
    cross = CrackCrossOneSegment(scrk1,xp1,xp2,norm,cross);
    return CrackCrossOneSegment(scrk2,xp1,xp2,norm,cross);
}

// Deepest Method to process one line segment in crack from particle scrk1 to scrk1->nextSeg
// Checks extents of that segment first. If fails, finally do line segment crossing algorithm
short CrackHeader::CrackCrossOneSegment(CrackSegment *scrk1,Vector *xp1,Vector *xp2,Vector *norm,short cross) const
{
    // next segment, exit if none (i.e., terminal particle)
    CrackSegment *scrk2 = scrk1->nextSeg;
    if(scrk2 == NULL) return cross;
    
    // check this segment's extents
	if(!LineIsInExtents(xp1,xp2,scrk1->cnear,scrk1->cfar)) return cross;
    
    // Now must check for crossing
	Vector xp3,xp4;
    double area1,area2;
    
    // first point, but if first point and it is exterior, extrapolate
	// backwards 4 times the segment length
    xp3.x=scrk1->cp.x;
    xp3.y=scrk1->cp.y;
	if(scrk1==firstSeg && scrk1->tipMatnum==EXTERIOR_CRACK)
    {   xp3.x-=4.*(scrk2->cp.x-xp3.x);
		xp3.y-=4.*(scrk2->cp.y-xp3.y);
	}
    
    // second point, but if last point and it is exterior, extrapolate
	// forward 4 times the segment length
    xp4.x=scrk2->cp.x;
    xp4.y=scrk2->cp.y;
    if(scrk2==lastSeg)
    {	if(scrk2->tipMatnum==EXTERIOR_CRACK)
        {	xp4.x-=4.*(xp3.x-xp4.x);
            xp4.y-=4.*(xp3.y-xp4.y);
        }
    }
        
	// checking areas of various triangles
	// See JAN0048-7 for details
    while(TRUE)
    {   // first two areas (123 and 124)
        area1=Triangle(xp1,xp2,&xp3);
        area2=Triangle(xp1,xp2,&xp4);
        
        // first area negative
        if(area1<0.)
        {	if(area2>0.)
            {   if(Triangle(&xp3,&xp4,xp1)<=0.) break;
            
                // TRUE mean - + + (- or 0) (0 means node on crack)
                if(Triangle(&xp3,&xp4,xp2)<=0.) goto above;
            }
            
            else if(area2==0.)
            {   if(Triangle(&xp3,&xp4,xp1)<=0.) break;
            
                // TRUE means - 0 + 0 (node on pt 4) or - 0 + - (pt 4 between mpt and node)
                if(Triangle(&xp3,&xp4,xp2)<=0.) goto above;
            }
        }
        
        // first area positive
        else if(area1>0.)
        {	if(area2<0.)
            {   if(Triangle(&xp3,&xp4,xp1)>=0.) break;
                
                // TRUE means + - - (+ or 0) (0 means node on crack)
                if(Triangle(&xp3,&xp4,xp2)>=0.) goto below;
            }
                
            else if(area2==0.)
            {   if(Triangle(&xp3,&xp4,xp1)>=0.) break;
                
                // TRUE means + 0 - 0 (node on pt 4) or + 0 - + (pt 4 between mpt and node)
                if(Triangle(&xp3,&xp4,xp2)>=0.) goto below;
            }
        }
        
        // first area zero
        else
        {	if(area2<0.)
            {   if(Triangle(&xp3,&xp4,xp1)>=0.) break;
                
                // TRUE means 0 - - 0 (node on pt 3) or 0 - - + (pt 3 between mpt and node) 
                if(Triangle(&xp3,&xp4,xp2)>=0.) goto below;
            }
                
            else if(area2>0.)
            {   if(Triangle(&xp3,&xp4,xp1)<=0.) break;
                
                // TRUE means 0 + + 0 (node on pt 3) or 0 + + - (pt 3 between mpt and node)
                if(Triangle(&xp3,&xp4,xp2)<=0.) goto above;
            }
        }
        
        // it does not cross
        break;
        
        // toggle the setting in case there are multiple crossings
    above:
        if(cross==NO_CRACK)
        {	cross=ABOVE_CRACK;
            norm->y=xp3.x-xp4.x;			// -¶x
            norm->x=xp4.y-xp3.y;			// ¶y
        }
        else
            cross=NO_CRACK;
        break;
    below:
        if(cross==NO_CRACK)
        {	cross=BELOW_CRACK;
            norm->y=xp3.x-xp4.x;			// -¶x
            norm->x=xp4.y-xp3.y;			// ¶y
        }
        else
            cross=NO_CRACK;
        break;
    }
	
	// needed for 3D crack coding
	norm->z = 0.;
    return cross;
}

// Determine if line from particle (x1,y1) to node (x2,y2) crosses this crack
// Return ABOVE_CRACK (1), BELOW_CRACK (2), or NO_CRACK (0) and crack normal in norm
// This method uses hierarchical crack in a binary tree
short CrackHeader::CrackCrossOnce(double x1,double y1,double x2,double y2,CrackSegment **crossSeg) const
{
    // recursive method to traverse tree hierarchy
	Vector xp1 = MakeVector(x1,y1,0.);
	Vector xp2 = MakeVector(x2,y2,0.);
    return CrackCrossLeafOnce(rootLeaf,&xp1,&xp2,crossSeg);
}

// Recursive Method to process each leaf in hierarchical traversal
short CrackHeader::CrackCrossLeafOnce(CrackLeaf *leaf,Vector *xp1,Vector *xp2,CrackSegment **crossSeg) const
{
    // check extents this leaf, return current cross if not in this leaf's extent
	if(!LineIsInExtents(xp1,xp2,leaf->cnear,leaf->cfar)) return NO_CRACK;
    
    // It is in extent of this leaf
    // if not terminal, go on to the child leaves
	short cross = NO_CRACK;
    if(!leaf->ChildrenAreSegments())
    {   CrackLeaf *child1,*child2;
        leaf->GetChildLeaves(&child1,&child2);
        cross = CrackCrossLeafOnce(child1,xp1,xp2,crossSeg);
		if(cross!=NO_CRACK) return cross;
        if(child2!=NULL) cross = CrackCrossLeafOnce(child2,xp1,xp2,crossSeg);
        return cross;
    }
    
    // Method 1: This code checks each segment now in a subroutine
	Vector norm;
    CrackSegment *scrk1,*scrk2;
    leaf->GetChildSegments(&scrk1,&scrk2);
    cross = CrackCrossOneSegment(scrk1,xp1,xp2,&norm,cross);
	if(cross!=NO_CRACK)
	{	*crossSeg = scrk1;
		return cross;
	}
    cross = CrackCrossOneSegment(scrk2,xp1,xp2,&norm,cross);
	if(cross!=NO_CRACK)
	{	*crossSeg = scrk2;
		return cross;
	}
	return NO_CRACK;
}

// When crack is first created at start of calculations, create all the CrackLeaf
// objects needed to describe the crack as a binary tree starting from
// the rootleaf
//
// Strategy
// 1. Create leafs each containing two neighboring segments (2D). In (3D) presort
//    to get segments in some proximity order and combine pairwise. Will have
//    number of segments/2 leaves (round up if last leaf has only one segment)
//    The defined crack must have at least two segments (double check this and cause error elsewhere if < 2)
// 2. Make parents of leaves by combining two leaves are current level. With
//    have number of leaves/2 parents (round up if last parent has only
//    one leaf.
// 3. When last level ended with only one parent, then done and that
//    parent becomes the root leaf in the crack tree
//
// Notes:
// 1. All cracks require 1 segment (or 2 points) to create a hierarchy. The code causes a fatal
//	  error if try to create crack with on one crack particle. The error occurs because this
//    method returns false
//
// throws std::bad_alloc
// CrackHeader3D overrides this
bool CrackHeader::CreateHierarchy(void)
{
    // all cracks have at least 2 points (and 1 segment)
    CrackSegment *scrk1 = firstSeg;
    if(scrk1==NULL) return false;
    CrackSegment *scrk2 = scrk1->nextSeg;
    if(scrk2==NULL) return false;
    CrackLeaf *firstLeaf=NULL,*lastLeaf=NULL,*leaf;
    int numLeaves=0;
    
    // This code converts ordered list of crack segments into a B-tree
    
    // Create terminal leafs each by pairwise groupin on crack particles or
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
	
    // Recursively create parent leaves by pairwise grouping leafs
    // at current level until the parent leval has only a single leaf
    while(numLeaves>1)
    {   // Pairwise grouping of current level that starts with firstLeaf
        // Note that numLeaves is number of leave in prior level and if that is
        //      1, we are done (ir it is the root level
        CrackLeaf *leaf1 = firstLeaf;
        CrackLeaf *leaf2 = leaf1->nextLeaf;
        firstLeaf = NULL;
        numLeaves = 0;
        
        // This creates the new level and counts number of leaves in that level
        // If only 1 leaf, it it root level
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
    return true;
}

// When crack moves, go through the hierarchy and recalculate the extents
// of each leaf in the crack.
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
// throws std::bad_alloc
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

// Determine if contour segment after nextPt crosses a segment of this crack
// Stops when finds first crossing. This would miss unlikely situation
//    where crack crosses the same contour segment more than once
CrackSegment *CrackHeader::ContourCrossCrack(ContourPoint *nextPt,Vector *crossPt) const
{
	// get line segment in the contour
    double x1,x2,y1,y2;
    switch(nextPt->orient)
    {	case HORIZONTAL:
            x1=nextPt->node->x;
            x2=nextPt->nextPoint->node->x;
            y1=nextPt->node->y;
            y2=y1;
			break;
		case VERTICAL:
            y1=nextPt->node->y;
            y2=nextPt->nextPoint->node->y;
            x1=nextPt->node->x;
            x2=x1;
			break;
		default:
			// error if grid not along x and y axes
			throw "The J Contour is not along x and y axes.";
	}
    
    // recursive call through crack hierarchy
    return ContourCrossLeaf(rootLeaf,x1,y1,x2,y2,crossPt,nextPt->orient);
}

// Recursive Method to process each leaf in hierarchical traversal
CrackSegment *CrackHeader::ContourCrossLeaf(CrackLeaf *leaf,double x1,double y1,double x2,double y2,Vector *crossPt,int orient) const
{
    // check extents this leaf, return NULL if not
	if(!LineIsInExtents(x1,y1,x2,y2,leaf->cnear,leaf->cfar)) return NULL;
    
    // It is in extent of this leaf
    
    // if not terminal leaf, go on to the child leaves
    CrackSegment *startSeg;
    if(!leaf->ChildrenAreSegments())
    {   CrackLeaf *child1,*child2;
        leaf->GetChildLeaves(&child1,&child2);
        startSeg = ContourCrossLeaf(child1,x1,y1,x2,y2,crossPt,orient);
        if(startSeg!=NULL) return startSeg;
        if(child2!=NULL) return ContourCrossLeaf(child2,x1,y1,x2,y2,crossPt,orient);
        return NULL;
    }
    
    // This code checks each segment in this terminal leaf
    CrackSegment *scrk1,*scrk2;
    leaf->GetChildSegments(&scrk1,&scrk2);
    if(SegmentsCross(scrk1,x1,y1,x2,y2,crossPt,orient)) return scrk1;
    if(SegmentsCross(scrk2,x1,y1,x2,y2,crossPt,orient)) return scrk2;
    return NULL;
}

// Determine if two line-segments cross
bool CrackHeader::SegmentsCross(CrackSegment *scrk1,double x1,double y1,double x2,double y2,Vector *crossPt,int orient) const
{
    // next segment, exit if none (i.e., terminal particle)
    CrackSegment *scrk2 = scrk1->nextSeg;
    if(scrk2 == NULL) return false;
    
    // check extents this segment, return NULL if not
	if(!LineIsInExtents(x1,y1,x2,y2,scrk1->cnear,scrk1->cfar)) return false;

    // find intersection and see if in the contour
    double dx,dy;
    if(orient==HORIZONTAL)
    {   dy = scrk2->cp.y-scrk1->cp.y;
        if(DbleEqual(dy,0.)) return (x1==scrk1->cp.x);		// parallel lines
        dx = scrk2->cp.x-scrk1->cp.x;
        crossPt->y = y1;
        crossPt->x = (crossPt->y-scrk1->cp.y)*dx/dy + scrk1->cp.x;
        if(x1<x2)
            return (crossPt->x>x1 && crossPt->x<=x2);
        else
            return (crossPt->x>x2 && crossPt->x<=x1);
    }
    
    // rest is VERTICAL
    dx = scrk2->cp.x-scrk1->cp.x;
    if(DbleEqual(dx,0.)) return (y1==scrk1->cp.y);         // parallel lines
    dy = scrk2->cp.y-scrk1->cp.y;
    crossPt->x = x1;
    crossPt->y = (crossPt->x-scrk1->cp.x)*dy/dx + scrk1->cp.y;
    if(y1<y2)
        return (crossPt->y>y1 && crossPt->y<=y2);
    else
        return (crossPt->y>y2 && crossPt->y<=y1);
}

// If current J contour crossed any cracks check those cracks for crossing
// and if it crosses, stop the grow at that crack
// Return relative change made in grow
double CrackHeader::AdjustGrowForCrossing(Vector *grow,CrackSegment *crkTip,double cSize,Vector *tipDir)
{
	// exit if not in the contour
	if(crossedCracks==NULL) return 1.0;
    if(crossedCracks->firstObject==NULL) return 1.0;
	
	// get growing line segment
	double x1 = crkTip->cp.x;
	double y1 = crkTip->cp.y;
	
	// add 1/4 cell
	double x2 = x1 + grow->x + 0.25*cSize*tipDir->x;
	double y2 = y1 + grow->y + 0.25*cSize*tipDir->y;
	
	// loop over crossed cracks
	double p = 1.0;
	CrackSegment *crossSeg;
	CrossedCrack *nextCross = (CrossedCrack *)crossedCracks->firstObject;
	while(nextCross!=NULL)
	{	CrackHeader *nextCrack = nextCross->crack;
		if(nextCrack->CrackCrossOnce(x1,y1,x2,y2,&crossSeg)!=NO_CRACK)
		{	// propagation path is x = x1 + p*grow->x, y = y1 + p*grow->y
			// line along the crack segment is x = crossSeg->x + s*dxs, y = crossSeg->y + s*dys
			CrackSegment *nextSeg = crossSeg->nextSeg;
			double dxs = nextSeg->cp.x-crossSeg->cp.x;
			double dys = nextSeg->cp.y-crossSeg->cp.y;
			
			// equate x and solve for s : s = (x1 - crossSeg->x + p*grow->x)/dxs
			//    But, if dxs==0 then x1 + p*grow->x = crossSeg->x to give p = (crossSeg->x-x1)/grow->x
			//		   if grow->x==0 too, parallel lines so do nothing
			// substitute s into equated y's and solve for p
			//	       p*(grow->y*dxs - grow->x*dys) = (crossSeg->y - y1)*dxs + (x1 - crossSeg->x)*dys
			//         But, if (grow->y*dxs - grow->x*dys)==0, parallel lines so do nothing
			//    But, if dxs=0 p already found above
			double cp = grow->y*dxs - grow->x*dys;			// will be zero if lines are parallel
			if(cp!=0.)
			{	if(dxs==0.)
					p = (crossSeg->cp.x-x1)/grow->x;
				else
					p = ((crossSeg->cp.y - y1)*dxs + (x1 - crossSeg->cp.x)*dys)/cp;
				
				// check last segment length
				CrackSegment *prevSeg = crkTip->nextSeg;
				if(prevSeg==NULL) prevSeg = crkTip->prevSeg;
				double pxs = crkTip->cp.x-prevSeg->cp.x;
				double pys = crkTip->cp.y-prevSeg->cp.y;
				double plen = sqrt(pxs*pxs+pys*pys);					// length of crack tip segment
				double glen = sqrt(grow->x*grow->x+grow->y*grow->y);	// length to proposed new tip segment
				if(PropagateTask::cellsPerPropagationStep>.7)
				{	int numSegs= (int)(2*(PropagateTask::cellsPerPropagationStep+.25));
					glen /= (double)numSegs;
				}
				double pmin = plen/glen<0.5 ? 0.5 : 0.0;
				
				// Adjust grow, but if previouis segment was not tool small
				if(p>=pmin)
				{	grow->x *= p;
					grow->y *= p;
					cout << "# Crack " << GetNumber() << " intersected crack " << nextCrack->GetNumber();
					cout << " at (" << x1 + grow->x << "," << y1+grow->y << ") fraction = " << p << endl;
				}
				else
					p = 1.0;
			}
			break;
		}
		nextCross = (CrossedCrack *)nextCross->GetNextObject();
	}
	return p;
}

#pragma mark CrackHeader: 2D ACCESSORS

// Find length of current crack
double CrackHeader::Length(void) const
{
    CrackSegment *scrk=firstSeg;
    double length=0.,lastx,lasty,dist2,x,y;
    
    // first point
    if(scrk==NULL) return length;
    lastx=(scrk->surf[0].x+scrk->surf[1].x)/2.;
    lasty=(scrk->surf[0].y+scrk->surf[1].y)/2.;
    scrk=scrk->nextSeg;
    
    // loop over crack points
    while(scrk!=NULL)
    {	// add segment length
    	x=(scrk->surf[0].x+scrk->surf[1].x)/2.;
        y=(scrk->surf[0].y+scrk->surf[1].y)/2.;
        dist2=(x-lastx)*(x-lastx) + (y-lasty)*(y-lasty);
        length+=sqrt(dist2);
        lastx=x;
        lasty=y;
        
        // next segment
        scrk=scrk->nextSeg;
    }
    return length;
}

// Find length of all segments that do not have a traction law
double CrackHeader::DebondedLength(void) const
{
    CrackSegment *scrk=firstSeg;
    double length=0.,lastx,lasty,dist2,x,y;
    
    // first point
    if(scrk==NULL) return length;
    lastx=(scrk->surf[0].x+scrk->surf[1].x)/2.;
    lasty=(scrk->surf[0].y+scrk->surf[1].y)/2.;
    scrk=scrk->nextSeg;
    
    // loop over crack points
    while(scrk!=NULL)
    {    // add segment length
        x=(scrk->surf[0].x+scrk->surf[1].x)/2.;
        y=(scrk->surf[0].y+scrk->surf[1].y)/2.;
        if(scrk->MatID()<0)
        {   dist2=(x-lastx)*(x-lastx) + (y-lasty)*(y-lasty);
            length+=sqrt(dist2);
        }
        lastx=x;
        lasty=y;
        
        // next segment
        scrk=scrk->nextSeg;
    }
    return length;
}

// given whichTip find crack tip segment
CrackSegment *CrackHeader::GetCrackTip(int whichTip) const
{	return (whichTip==START_OF_CRACK) ? firstSeg : lastSeg ;
}

// given whichTip find crack tip segment
CrackSegment *CrackHeader::GetAdjToCrackTip(int whichTip) const
{	return (whichTip==START_OF_CRACK) ? firstSeg->nextSeg : lastSeg->prevSeg ;
}

// given crack tip segment find which tip
int CrackHeader::GetWhichTip(CrackSegment *crkTip) const
{	return (crkTip==firstSeg) ? START_OF_CRACK  : END_OF_CRACK ;
}

// Deterimine what needs this crack tip has to do propagation calculation
// Return 0 (no need), NEED_JANDK, or NEED_J
int CrackHeader::CriterionNeeds(void) const
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
bool CrackHeader::GetAllowAlternate(int crkTipIdx) const { return allowAlternate[crkTipIdx]; }
void CrackHeader::SetAllowAlternate(int crkTipIdx,bool setting) { allowAlternate[crkTipIdx]=setting; }

// describe velocity field
void CrackHeader::Describe(void) const
{	cout << "# crack#=" << number << " thickness=" << thickness << " custom contact mat="
			<< customContactLawID << " has traction=" << hasTractionLaws << endl;
	cout << "#    from (" << firstSeg->cp.x << "," << firstSeg->cp.y << ") to ("
			<< lastSeg->cp.x << "," << lastSeg->cp.y << "), length = " << Length() ;
	cout << ", segments = " << NumberOfSegments() << endl;
}

// crack thickness for planar problems
void CrackHeader::SetThickness(double thick) { thickness = thick; }
double CrackHeader::GetThickness(void) const { return thickness; }
double *CrackHeader::GetThicknessPtr(void) { return &thickness; }

#pragma mark CrackHeader: ACCESSORS

// Set segments to number of keypoints to current number of segments
// Only used by 3D cracks
void CrackHeader::SetNumberKeypoints(void) { numberKeypoints = numberSegments; }

// Count Segments
int CrackHeader::NumberOfSegments(void) const { return numberSegments; }

// Count Facets (on for 3D cracks)
int CrackHeader::NumberOfFacets(void) const { return numberFacets; }

// crack number
void CrackHeader::SetNumber(int newNum) { number=newNum; }
int CrackHeader::GetNumber(void) const { return number; }

// make it a fixed crack
void CrackHeader::SetFixedCrack(int fixSetting) { fixedCrack=fixSetting; }

// set custom contact law ID
void CrackHeader::SetContactLawID(int newID) { customContactLawID = newID; }

// set/get custom traction propagate low ID (1 based)
void CrackHeader::SetTractionPropID(int newID) { customTractionPropID = newID; }
int CrackHeader::GetTractionPropID(void) { return customTractionPropID; }

// return if any segment has traction law material
bool CrackHeader::GetHasTractionLaws(void) const { return hasTractionLaws; }

// Type of crack
bool CrackHeader::IsThreeD(void) const { return numberFacets>0; }

#pragma mark CrackHeader: Class methods

// Signed area of a 2D triangle from coordinates
double CrackHeader::Triangle(double x1,double y1,double x2,double y2,double x3,double y3)
{   return(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));
}

// Signed area of a triangle in x-y plane from vectors
double CrackHeader::Triangle(Vector *xp1,Vector *xp2,Vector *xp3)
{   return(xp1->x*(xp2->y-xp3->y) + xp2->x*(xp3->y-xp1->y) + xp3->x*(xp1->y-xp2->y));
}

// Find location for spline interpolation on crack surfaces
// We fit the first four crack particle to a B spline. The value of t (from 0 to 4)
//   is location relative to crack tip used to find crack opening displacement. The spline
//   curve has four curves defined by these particles. The default method is to use
//   the crack particles as the control points. An option called CUBIC_INTERPOLATION changes
//   crack particle to interpolation points (and control points are found. An options
//   called LINEAR_INTERPOLATION does not use splines (or setting below
// The B spline curve used depends on t (0 to 1 for first, to 2 for second, to 3 for third)
void CrackHeader::SetCodLocation(double t)
{
	// we do not allow t outside 0 to 3
	if(t<0.) t=0.;
	if(t>3.) t=3.;
	
	// Pick one of first three B spline curves
	if(t<=1.)
		codInterval=0;
	else if(t<=2.)
	{	codInterval=1;
		t-=1.;				// range 0 to 1 in the interval
	}
	else
	{	codInterval=2;
		t-=2.;				// range 0 to 1 in the interval
	}
	codLocation=t;			// location 0 to 1 in the choosen curvne
	
	// coefficients to P(t) curve
	bezArg[0]=(1.-t)*(1.-t)*(1.-t);
	bezArg[1]=3.*t*(1.-t)*(1.-t);
	bezArg[2]=3.*t*t*(1.-t);
	bezArg[3]=t*t*t;
	
	// coefficient for P'(t)
	// These could be used to find crack tangent and normal, but current is not used.
	bezDer[0]=-3.*(1.-t)*(1.-t);
	bezDer[1]=3.*(1.-t)*(1.-t) - 6.*t*(1.-t);
	bezDer[2]=6.*t*(1-t) - 3.*t*t;
	bezDer[3]=3.*t*t;
	
	// current code settings has t=1 on input or
	// codLocation = 1, codInterval=0
	// bezArg = {0,0,0,1}
	// bezDer = {0,0,-3,3}
}

// Check is lines is within extent of leaf or segment
// return false if line within extents or false otherwise
// See JANOSU-6-66
bool CrackHeader::LineIsInExtents(Vector *xp1,Vector *xp2,double *cnear,double *cfar)
{	return LineIsInExtents(xp1->x,xp1->y,xp2->x,xp2->y,cnear,cfar);
}
bool CrackHeader::LineIsInExtents(double x1,double y1,double x2,double y2,double *cnear,double *cfar)
{
	// Pi = (1,0)
    if(fmax(x1,x2) < cnear[0]) return false;
    if(fmin(x1,x2) > cfar[0]) return false;

	// Pi = (0,1)
    if(fmax(y1,y2) < cnear[1]) return false;
    if(fmin(y1,y2) > cfar[1]) return false;
    
	// Pi = (1,1)
    double Pib = x1 + y1;                       // Pi.b
    double Pia = x2 + y2 - Pib;                 // Pi.a with Pi = (1,1)
    if(Pia>0.)
    {   if(cnear[2]-Pib > Pia) return false;
        if(cfar[2] < Pib) return false;
    }
    else
    {   // This works for Pia=0 as well
        if(cfar[2]-Pib < Pia) return false;
        if(cnear[2] > Pib) return false;
    }
    
	// Pi = (1,-1)
    Pib = x1 - y1;                              // Pi.b
    Pia = x2 - y2 - Pib;                        // Pi.a with Pi = (1,-1)
    if(Pia>0.)
    {   if(cnear[3]-Pib > Pia) return false;
        if(cfar[3] < Pib) return false;
    }
    else
    {   // This works for Pia=0 as well
        if(cfar[3]-Pib < Pia) return false;
        if(cnear[3] > Pib) return false;
    }
	
	// passes all tests
	return true;
}

// Function to prehash cracks for faster crossing algorithm
#ifdef PREHASH_CRACKS
void CrackHeader::UpdateElementCrackList(int cracknum)
{
	// need number of elements 
	int horiz, vert, depth;
	mpmgrid.GetGridPoints(&horiz, &vert, &depth);
	
	// PHCNew
	// GetGridPoints() get number of nodes, subtract one
	// to get number of elements in each row (horiz) and number of rows (vert)
	horiz--; vert--;

	// Get starting segments
	CrackSegment *CrackPnts = firstSeg;
    
	// PHCNew
	int Element2,Element1 = CrackPnts->planeElemID(); //zero based
	int y2,y1 = Element1 / horiz;			// row 0 to vert-1
	int x2,x1 = Element1%horiz;				// col 0 to horiz-1

	// loop over crack segments
	while (CrackPnts->nextSeg != NULL) {

		// PHCNew
		Element2 = CrackPnts->nextSeg->planeElemID();		// element of second point

		// get coordinates +1 for GIMP and CPDI
		y2 = Element2 / horiz;					// row 0 to vert-1
		x2 = Element2%horiz;					// col 0 to horiz-1
		int maxy, maxx, miny, minx;

		// get bounding box
		// get max/min of y plus 1;
		if (y1 > y2) {
			maxy = y1 + 1;
			miny = y2 - 1;
		}
		else {
			maxy = y2 + 1;
			miny = y1 - 1;
		}

		// get max/min of x plus 1;
		if (x1 > x2) {
			maxx = x1 + 1;
			minx = x2 - 1;
		}
		else {
			maxx = x2 + 1;
			minx = x1 - 1;
		}

		// make sure we are not beyond edges
		// PHCNew - do need to stop 1 below horiz
		maxx = fmin(maxx, horiz - 1);
		maxy = fmin(maxy, vert - 1);
		minx = fmax(minx, 0);
		miny = fmax(miny, 0);
		
		// loop over elements in the bounding box and add this crack to their seen list
		for (int j = miny; j <= maxy; j++) {
			int rowStart = j*horiz;
			for (int i = minx; i <= maxx; i++) {
				theElements[rowStart + i]->PushCrackNumOnList(cracknum);
			}
		}

		// PHCNew
		// next segment
		Element1 = Element2;
		x1 = x2;
		y1 = y2;
		CrackPnts = CrackPnts->nextSeg;
	}
}
#endif
