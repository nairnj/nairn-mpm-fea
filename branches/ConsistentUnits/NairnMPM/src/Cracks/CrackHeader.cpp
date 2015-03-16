/********************************************************************************
    CrackHeader.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Apr 3 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

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
#include "Nodes/NodalPoint2D.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Cracks/CrackLeaf.hpp"
#include "Cracks/CrossedCrack.hpp"
#include "Read_XML/ParseController.hpp"
#include "Custom_Tasks/PropagateTask.hpp"

// Include to store J using 1 term in J2. Calculation should use
// two terms to get that result in J1
//#define RECORD_ONE_AND_TWO_TERM_RESULTS

using namespace std; 

// class statics
double CrackHeader::codLocation;
int CrackHeader::codInterval;
double CrackHeader::bezArg[4];
double CrackHeader::bezDer[4];
int CrackHeader::warnNodeOnCrack;
int CrackHeader::warnThreeCracks;

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
    firstSeg=NULL;
    lastSeg=NULL;
    numberSegments=0;
	fixedCrack=FALSE;
	customContact=FALSE;
	hasTractionLaws=FALSE;
	thickness=1.0;				// for crack tip heating and tractions in mm, will default to grid thickness if set
	allowAlternate[0]=allowAlternate[1]=TRUE;
	crossedCracks = NULL;
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
    }
    else
	{	// no need to add a zero length segment
		if(DbleEqual(cs->x,lastSeg->x) && DbleEqual(cs->y,lastSeg->y))
		{	// but maybe want new tip material or traction law
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
    4. Transfer Crack tip properties to new crack tip
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
    
    ExtendHierarchy(cs);
    
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
		double fn[maxShapeNodes],shapeNorm;
		int j,nodeCounter;
		Vector delv,cpos,vcm;
		int nds[maxShapeNodes],numnds;
		
		// loop over crack points
		while(scrk != NULL)
		{	if(!fixedCrack)
			{	// get element and shape functinos
				iel=scrk->planeInElem-1;			// now zero based
				cpos.x=scrk->x;
				cpos.y=scrk->y;
				theElements[iel]->GetShapeFunctionsForCracks(&numnds,fn,nds,&cpos);
				
				// initialize
				ZeroVector(&delv);
				nodeCounter=0;
				shapeNorm=0.;
				
				// extrapolate to particle
				for(j=1;j<=numnds;j++)
				{	if(nd[nds[j]]->GetCMVelocityTask8(&vcm))
					{	AddScaledVector(&delv,&vcm,fn[j]);
						nodeCounter++;
						shapeNorm+=fn[j];
					}
				}
				
				// move it or collapse it
				if(nodeCounter>0)
				{	// renormalize and multiply by dt to get displacement
					ScaleVector(&delv,timestep/shapeNorm);
					scrk->MovePosition(delv.x,delv.y);		// in mm
					
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
					
					// development flag to collapse wide open cracks during cutting
					if(fmobj->dflag[0]==4 && scrk->MatID())
					{	// find COD
						double codx=scrk->surfx[0]-scrk->surfx[1];
						double cody=scrk->surfy[0]-scrk->surfy[1];
						double cod=sqrt(codx*codx+cody*cody);
						if(cod>0.75)
						{	scrk->x=(scrk->surfx[0]+scrk->surfx[1])/2.;
							scrk->y=(scrk->surfy[0]+scrk->surfy[1])/2.;
							if(!scrk->FindElement()) return false;
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
    double fn[maxShapeNodes],surfaceMass;
    short js=side-1,nodeCounter,j;
	int numnds,nds[maxShapeNodes];
    Vector delv,cpos;
    
    // loop over crack points
    while(scrk!=NULL)
	{	if(!fixedCrack)
		{	// get element
			iel = scrk->surfInElem[js]-1;			// now zero based
			cpos.x = scrk->surfx[js];
			cpos.y = scrk->surfy[js];
			theElements[iel]->GetShapeFunctionsForCracks(&numnds,fn,nds,&cpos);
            
			// initialize
			ZeroVector(&delv);
			surfaceMass = 0;
			nodeCounter = 0;
			
			// extrapolate those with velocity to the particle
			for(j=1;j<=numnds;j++)
			{	if(nd[nds[j]]->IncrementDelvSideTask8(side,number,fn[j],&delv,&surfaceMass,scrk))
					nodeCounter++;
			}

			// if CRACK_SURFACE_BY_MOMENTUM_EXTRAP is defined
			//     delv is Sum(fi pi) and surfaceMass = Sum(fi mi)
			// otherwise
			//     elv is Sum(fi vi) = Sum(fi pi/mi) and surfaceMass = Sum(fi)
			// Both normalize to get velocity and multiply by dt to get displacement
			if(nodeCounter>0) ScaleVector(&delv,timestep/surfaceMass);
			
			// this method does not normalize shape functions
			//ScaleVector(&delv,timestep);
           
			// move it (if returns true, check location of other side for element move again because it moved too)
			if(scrk->MoveSurfacePosition(side,delv.x,delv.y,(nodeCounter>0)))		// in mm
			{	if(!scrk->FindElement(ABOVE_CRACK)) return false;
			}
			
			// did surface move elements
			if(!scrk->FindElement(side)) return false;
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
{
	if(!hasTractionLaws) return;
	CrackSegment *cs=firstSeg;
	while(cs!=NULL)
	{	cs->AddTractionForceSeg(this);
		cs=cs->nextSeg;
	}
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
    // This is method that is currently active
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
	Vector crackDir,crossPt;
	bool secondTry;
	
    /* Calculate J-integrals for the ith crack tip
    */
	
	// Clear previously cross cracks list
	// TNote: this linked list is for cracks crossed by the contours. It has cracks
	// from both tips in one list. The propagation task can look at these cracks to
	// check for propagation across another crack. The fact that both tips are in one
	// list is not an issue, because the list will be short and checking will be
	// efficient even if a lot of cracks in the simulation
	if(crossedCracks!=NULL)
	{	// check if this crack is already in the list
		CrossedCrack *nextCross = (CrossedCrack *)crossedCracks->firstObject;
		while(nextCross!=NULL)
		{	CrossedCrack *hold = (CrossedCrack *)nextCross->GetNextObject();
			delete nextCross;
			nextCross = hold;
		}
		delete crossedCracks;
		crossedCracks = NULL;
	}
	
	// it may try two contours at each crack tip. First try is at NearestNode().
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
		crackPt=NULL;					// first crack pt in the contour
		try
		{
			/* Task 2: find ccw nodal points JGridSize from crack tip nodal point
				Find orientation of each line segment
			*/
			gridElem=tipCrk->planeInElem-1;
			gridNode=theElements[gridElem]->NearestNode(tipCrk->x,tipCrk->y,&nextNearest);
			if(secondTry) gridNode=nextNearest;
			
			// set material type based on material near the tip.
			int oldnum = tipCrk->tipMatnum;
			tipCrk->FindCrackTipMaterial(oldnum);
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
            //        2
            //    ---------
            //   |         |
            // 3 |    *    | 1  * is the crack tip
            //   |         |
            //    ----|----
            //    4     0
			prevPt=crackPt=new ContourPoint(nd[gridNode]);
			int numPts=JGridSize;
			double cxmin=9e99,cxmax=-9e99,cymin=9e99,cymax=-9e99;
			for(j=0;j<5;j++)
			{   gridNode=theElements[gridElem]->NextNode(gridNode);
				nextPt=new ContourPoint(nd[gridNode]);
				prevPt->SetNextPoint(nextPt);
				prevPt=nextPt;
				for(i=1;i<numPts;i++)
				{   gridElem=theElements[gridElem]->Neighbor(gridNode);
					if(gridElem<0)
						throw "J integral contour at a crack tip does not fit in the grid";
					gridNode=theElements[gridElem]->NextNode(gridNode);
					nextPt=new ContourPoint(nd[gridNode]);
					prevPt->SetNextPoint(nextPt);
					prevPt=nextPt;
				}
				
				// check corners for extent of contour rectangle
				cxmin=min(cxmin,nd[gridNode]->x);
				cxmax=max(cxmax,nd[gridNode]->x);
				cymin=min(cymin,nd[gridNode]->y);
				cymax=max(cymax,nd[gridNode]->y);

                // adjust to half edge when about to do j=4 edge
				numPts = (j==3) ? JGridSize-1 : 2*JGridSize;
			}
			// connect end to start
			prevPt->SetNextPoint(crackPt);
			
			// Task 3 was to find crossing point, now done below

 			/* Task 4: Loop over all segments and evaluate J integral (from the primary term)
				Transform to crack plane and save results
			*/
			DispField *sfld1,*sfld2;
			// save point, crack segement, and contour segment where the crack crosses the target crack
			CrackSegment *startSeg = NULL;
			ContourPoint *crossContourPt = NULL;
			bool crossesOtherCracks = false;

			// Initialize contour integration terms
			Jx1=Jy1=0.0;				// J-integral components from the first term
			double tractionEnergy=0.,bridgingReleased=0.;		// for traction laws
			nextPt=crackPt;				// start tracing the contour
			int count=0;				// particles in the nodal fields
			double r1 = 1.,r2 = 1.;		// for axisymmetric J's
			
			// loop over nodes in the contour
			while(true)
			{   // J integral along the countour line from node1 to node2
				NodalPoint *node1=nextPt->node;
				NodalPoint *node2=nextPt->nextPoint->node;

				// does line from node1 to node2 cross this crack&
				Vector crossPt1;
				int crossCount = 0;
				int crackNum = 0;
				double fract;
				
				// check for crack crossing, but if last segment crosses, the next pass
				// through the loop will start with a phantomNode and therefore do not
				// do the check. Note that do not need to check second node as phantom node
				// because that never occurs at this stage.
				if(!nextPt->phantomNode)
				{	// Check crossing of this crack first
					CrackSegment *foundSeg = ContourCrossCrack(nextPt,&crossPt1);
					if(foundSeg != NULL)
					{	crossCount++;
						if(crossCount==2)
						{   // error unless new crossPt is the same, which implies endpoints for two adjacent segments
							// two identical endpoints are accepted, but otherwise an error
							if(!(DbleEqual(crossPt.x,crossPt1.x) && DbleEqual(crossPt.y,crossPt1.y)))
							{	if(secondTry)
                                {   throw "Two different crossings between J-path and a crack";
                                }
                                else
                                    throw "";
                            }
						}
						else if(crossCount>2)
						{   // only gets here if found endpoints before and now cannot be another matching endpoint
							throw "Two different crossings between J-path and a crack";
						}
						else
						{   // save crossing point on target crack and segment
							crossContourPt = nextPt;
							crossPt.x = crossPt1.x;
							crossPt.y = crossPt1.y;
							// find crack particle closer to the crack tipstart
							startSeg = crkTipIdx==END_OF_CRACK ? foundSeg->nextSeg : foundSeg ;
							crackNum = GetNumber();
							fract=nextPt->Fraction(crossPt);
						}
					}
					
					// if does not cross the target crack, look for an interating crack
					// This stops at first intersecting crack and thus would miss more than
					//   one other crack intersecting the same segment.
					if(foundSeg == NULL)
					{	CrackHeader *nextCrack = firstCrack;
						while(nextCrack!=NULL)
						{	if(nextCrack != this)
                            {	foundSeg = nextCrack->ContourCrossCrack(nextPt,&crossPt1);
                                if(foundSeg != NULL)
                                {	// found crack crossing
                                    crackNum = nextCrack->GetNumber();
                                    crossesOtherCracks = true;
                                    fract=nextPt->Fraction(crossPt1);
									
									// keep list of crossed cracks
									CrossedCrack *nextCross = NULL;
									if(crossedCracks==NULL)
										crossedCracks = new ParseController();
									else
									{	// check if this crack is already in the list
										nextCross = (CrossedCrack *)crossedCracks->firstObject;
										while(nextCross!=NULL)
										{	if(nextCross->crack==nextCrack) break;
											nextCross = (CrossedCrack *)nextCross->GetNextObject();
										}
									}
									if(nextCross==NULL)
									{	crossedCracks->AddObject(new CrossedCrack(nextCrack));
									}
                                }
							}
                            nextCrack = (CrackHeader *)nextCrack->GetNextObject();
						}
					}
					
					// If this contour segment crosses a crack then must break into two contour segmnets
					// othersize use field [0] on both nodes
					if(crackNum!=0)
					{	// Create a phantom node at the crossing point
						NodalPoint *phantom=new NodalPoint2D(crackNum,crossPt1.x,crossPt1.y);
						phantom->PrepareForFields();
						
						// J integral will use phatom[0] up to a phantom contour point and phantom[1] starting from it by:
						//    Interpolate [0] from node1 and [1] or [2] that crosses crackNum from node2 to phantom [0]
						//    Interpolate [0] from node2 and [1] or [2] that crosses crackNum from node1 to phantom [1]
						phantom->Interpolate(node1,node2,fract,crackNum);
						node2 = phantom;
						
						// insert new ContourPoint with the phantom node
						ContourPoint *insertPt=new ContourPoint(phantom);
						insertPt->SetPhantomNode(true);
						insertPt->SetNextPoint(nextPt->nextPoint);
						nextPt->SetNextPoint(insertPt);
					}
				}
				
				// Get fields
				count += node1->GetFieldForCrack(nextPt->phantomNode,true,&sfld1,0);
				count += node2->GetFieldForCrack(nextPt->nextPoint->phantomNode,false,&sfld2,0);
#ifdef CONTOUR_PARTS
                cout << "#nodes " << node1->num << " to " << node2->num << ": ";
#endif

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
#ifdef CONTOUR_PARTS
                cout << "(Jxs,Jys)=(" << Jxs << "," << Jys << ")";
                cout << "(Jx1,Jy1)=(" << Jx1 << "," << Jy1 << ")" << endl;
#endif

				// on to next segment (switch field at mid point)
				nextPt=nextPt->nextPoint;
				if(nextPt==crackPt) break;
			}
			
			// if did not find any velocity fields in the countour, then the entire countour is
			// in empty space and J will be zero (and meaningless)
			if(count==0)
			{	throw "Section of the J Integral contour was in empty space";
			}

			// if startSeg is still NULL, then did not find any crossings
			if(startSeg == NULL)
			{	throw "A crack does not cross its J path";
			}
			//PrintContour(crackPt,crossContourPt,crossPt);
			
#ifdef PRINT_CROSS_STATUS
			if(crossesOtherCracks)
			{
                cout << "#... contour crosses other cracks" << endl;
			}
			else
            {
				cout << "#... contour crosses no other cracks" << endl;
            }
#endif

			/* Task 5: Evaluate J integral from the additional terms (GYJ)
				if requested */
			Jx2 = Jy2 = JxAS2 = 0.;
			if(JTerms==2)
			{   double rho,xp,yp,carea;
				double ax,ay,duxdx,duydx,duxdy,duydy;
				double vx,vy,dvxdx,dvydy,dvxdy,dvydx;
				double f2ForJx=0.,f2ForJy=0.,f2axisym=0.;
				count=0;	// number of particles within J-integral contour

				// integrate nonrigid particles
				for(int p=0;p<nmpmsNR;p++)
				{	xp=mpm[p]->pos.x;
					yp=mpm[p]->pos.y;
					if(xp>=cxmin && xp<cxmax && yp>=cymin && yp<cymax)
					{   // (xp,yp) in the contour
						count++;
						
						// Mass density g/mm^3
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
						
						// increment the integrands (g/mm^3)(mm/sec^2) = uN/mm^3 
						f2ForJx += rho*((ax*duxdx+ay*duydx)-(vx*dvxdx+vy*dvydx));
						f2ForJy += rho*((ax*duxdy+ay*duydy)-(vx*dvxdy+vy*dvydy));
						
						if(fmobj->IsAxisymmetric())
						{	// in axisymmetrix z is theta direction, etheta = u/r. but w=0, az=vz=0
							// Since w=0, no change to above terms, but have some static terms for Jx=Jr only
							// Stress Units uN/mm^2 (mm^3/g)
							Tensor sp = mpm[p]->ReadStressTensor();
							
							// Units uN/mm^3
							if(JContourType == AXISYM_BROBERG_J)
							{	// See Broberg, Cracks and Fraction (1999), page 65
								f2axisym += rho*(sp.xx*duxdx - sp.zz*mpm[p]->GetDwDz() + sp.xy*duydx)/xp;
							}
							else
							{	// Bergkvist and Huong called J3D/(a dphi)
								f2axisym += rho*(mpm[p]->GetWorkEnergy() - sp.zz*mpm[p]->GetDwDz())/crackr;
							}
						}
					}
				}
				
				if(count==0)
					throw "J Integral contour contains no particles";
				carea = (cxmax-cxmin)*(cymax-cymin)/count;		// area per particle in mm
				Jx2 = f2ForJx*carea;							// Jx2 in uN mm/mm^2 now
				Jy2 = f2ForJy*carea;							// Jy2 in uN mm/mm^2 now
				JxAS2 = f2axisym*carea;							// JxAS2 (for Jr in axisymmetric) in uN mm/mm^2
#ifdef CONTOUR_PARTS
				cout << "#...(Jx2,Jy2,JxAS2)=(" << Jx2 << "," << Jy2 << "," << JxAS2 << ")" << endl;
#endif
			}
			
			/* Task 6: Subtract energy due to tractions or for cracks in contact, subtract
				energy associated with shear stress (later not yet implemented thought)
			*/
			if(hasTractionLaws)
			{	CrackSegment *tipSegment = startSeg;
				tractionEnergy=startSeg->TractionEnergy(&crossPt,crkTipIdx,true,&tipSegment);
				bridgingReleased=startSeg->TractionEnergy(&crossPt,crkTipIdx,false,NULL);
				
				// extra traction correction if Bergkvist and Huong axisymmetric J integral
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
						previousPt.x = tipSegment->x;
						previousPt.y = tipSegment->y;
						closerSeg = (crkTipIdx==START_OF_CRACK) ? closerSeg->prevSeg : closerSeg->nextSeg ;
					}
					double extra = 0.;
					while(closerSeg!=NULL)
					{	double dr = closerSeg->x-previousPt.x;
						previousPt.x = closerSeg->x;
						previousPt.y = closerSeg->y;
						extra += closerSeg->TractionEnergy(&previousPt,crkTipIdx,true,NULL)*dr;
						closerSeg = (crkTipIdx==START_OF_CRACK) ? closerSeg->prevSeg : closerSeg->nextSeg ;
					}
					tractionEnergy += extra/crackr;
				}
			}
			else
			{	// set traction energy to energy due to shear if in contact (not implemented yet)
				tractionEnergy=0.;
				bridgingReleased=0.;
			}
			
			// add the two terms N mm/mm^2
			Jx = Jx1 + Jx2 - JxAS2;
			Jy = Jy1 + Jy2;
#ifdef CONTOUR_PARTS
			cout << "#...(Jx,Jy,Traction)=(" << Jx << "," << Jy << "," << tractionEnergy << ")";
			cout << "(cp.x,cp.y)=(" << crossPt.x << "," << crossPt.y << ")" << endl;
#endif
 
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
#ifdef RECORD_ONE_AND_TWO_TERM_RESULTS
			tipCrk->Jint.y = Jx1*crackDir.x + Jy1*crackDir.y - tractionEnergy;		// J by one term (temporary)
#endif
			tipCrk->Jint.z = tipCrk->Jint.x + bridgingReleased;						// Jrel or energy released in current state
#ifdef CONTOUR_PARTS
			cout << "#...J = " << tipCrk->Jint.x << endl;
#endif
			
			// end of try block on J calculation
			secondTry=FALSE;
		}
		catch(const char *msg)
		{	// throwing "" signals to try again with the next nearest node
			if(strlen(msg)==0)
			{	secondTry=TRUE;
                // contour is released below before trying again
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
		
		if(!secondTry) crkTipIdx++;
    } // end loop over crack tips
}

// print the contour (for debegging)
// Contour starts or crackPt, and crosses at crossPt in contour segement after prevPt
void CrackHeader::PrintContour(ContourPoint *crackPt,ContourPoint *crossContourPt,Vector &crossPt)
{
	cout << "# J Contour Nodes (node,orient,ds,nx,ny,phantom)" << endl;
	ContourPoint *nextPt = crackPt;
	while(true)
	{	cout << "#   " << nextPt->node->num << "," << nextPt->orient << ",";
		cout << nextPt->ds << "," << nextPt->norm.x << "," << nextPt->norm.y;
		cout << "," << nextPt->phantomNode;
		cout << endl;
		nextPt = nextPt->nextPoint;
		if(nextPt==crackPt) break;
	}
	if(crossContourPt!=NULL)
	{	cout << "#   Crossing Point: (" << crossPt.x << "," << crossPt.y << ")";
		cout << " after node " << crossContourPt->node->num << endl;
	}
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

// Add crack tip heating for all points in this crack
void CrackHeader::CrackTipHeating(void)
{
    CrackSegment *scrk=firstSeg;
	int iel;
	double fn[maxShapeNodes];
	Vector cpos;
	int numnds,i,nds[maxShapeNodes];
    
	// exit if no segments
    if(scrk==NULL) return;
    
    // loop over crack tips
    while(scrk!=NULL)
	{	iel=scrk->planeInElem-1;		// now zero based
		cpos.x=scrk->x;
		cpos.y=scrk->y;
		theElements[iel]->GetShapeFunctionsForCracks(&numnds,fn,nds,&cpos);
	
        // Add crack particle heating to each node in the element
        for(i=1;i<=numnds;i++)
		{	nd[nds[i]]->fcond+=scrk->HeatRate()*fn[i];
		}
		
		// next segment
        scrk=scrk->nextSeg;
	}
}

/* see if node is within tolerance (cell dimensions) of the crack tip
	needs regular grid with equal elements (could revise when needed), otherwise returns false
*/
bool CrackHeader::NodeNearTip(NodalPoint *ndi,double tol)
{
	if(!mpmgrid.IsStructuredGrid()) return FALSE;
    CrackSegment *scrk=firstSeg;
	double dx = (ndi->x-scrk->x)/mpmgrid.gridx;
	double dy = (ndi->y-scrk->y)/mpmgrid.gridy;
	double dist = sqrt(dx*dx+dy*dy);
	if(dist<tol) return true;
    scrk=lastSeg;
	dx = (ndi->x-scrk->x)/mpmgrid.gridx;
	dy = (ndi->y-scrk->y)/mpmgrid.gridy;
	dist = sqrt(dx*dx+dy*dy);
	if(dist<tol) return true;
	return false;
}

#pragma mark GLOBAL EXTENT CRACKS

// Determine if line from particle (1) to node (2) crosses this crack
// Return ABOVE_CRACK (1), BELOW_CRACK (2), or NO_CRACK (0) and crack normal in norm
// This method uses global extents for crack, if that fails it checks all
//      segments
// This method is no longer used, but is left here to check the hierarchical scheme
//      It is called from CFFlatCrossing()
short CrackHeader::FlatCrackCrossTest(double x1,double y1,double x2,double y2,Vector *norm)
{
    double x3,y3,x4,y4;
    CrackSegment *scrk=firstSeg;
    short cross=NO_CRACK;
    double area1,area2;
    
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
    // Points: 1 is particle, 2 is node, segment from 3 to 4
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

#pragma mark HIERARCHICAL CRACKS

// Determine if line from particle (x1,y1) to node (x2,y2) crosses this crack
// Return ABOVE_CRACK (1), BELOW_CRACK (2), or NO_CRACK (0) and crack normal in norm
// This method uses hierarchical crack in a binary tree
short CrackHeader::CrackCross(double x1,double y1,double x2,double y2,Vector *norm) const
{
    // recursive method to traverse tree hierarchy
    return CrackCrossLeaf(rootLeaf,x1,y1,x2,y2,norm,NO_CRACK);
}

// Recursive Method to process each leaf in hierarchical traversal
short CrackHeader::CrackCrossLeaf(CrackLeaf *leaf,double x1,double y1,double x2,double y2,Vector *norm,short cross) const
{
    // check extents this leaf, return current cross if not in this leaf's extent
	if(!LineIsInExtents(x1,y1,x2,y2,leaf->cnear,leaf->cfar)) return cross;
    
    // It is in extent of this leaf
    // if not terminal, go on to the child leaves
    if(!leaf->ChildrenAreSegments())
    {   CrackLeaf *child1,*child2;
        leaf->GetChildLeaves(&child1,&child2);
        cross = CrackCrossLeaf(child1,x1,y1,x2,y2,norm,cross);
        if(child2!=NULL) cross = CrackCrossLeaf(child2,x1,y1,x2,y2,norm,cross);
        return cross;
    }
    
    // Method 1: This code checks each segment now in a subroutine
    CrackSegment *scrk1,*scrk2;
    leaf->GetChildSegments(&scrk1,&scrk2);
    cross = CrackCrossOneSegment(scrk1,x1,y1,x2,y2,norm,cross);
    return CrackCrossOneSegment(scrk2,x1,y1,x2,y2,norm,cross);
}

// Deepest Method to process one line segment in crack from particle scrk1 to scrk1->nextSeg
// Checks extents of that segment first. If fails, finally do line segment crossing algorithm
short CrackHeader::CrackCrossOneSegment(CrackSegment *scrk1,double x1,double y1,double x2,double y2,Vector *norm,short cross) const
{
    // next segment, exit if none (i.e., terminal particle)
    CrackSegment *scrk2 = scrk1->nextSeg;
    if(scrk2 == NULL) return cross;
    
    // check this segment's extents
	if(!LineIsInExtents(x1,y1,x2,y2,scrk1->cnear,scrk1->cfar)) return cross;
    
    // Now must check for crossing
    double x3,y3,x4,y4;
    double area1,area2;
    
    // first point, but if first point and it is exterior, extrapolate
	// backwards 4 times the segment length
    x3=scrk1->x;
    y3=scrk1->y;
	if(scrk1==firstSeg && scrk1->tipMatnum==EXTERIOR_CRACK)
    {   x3-=4.*(scrk2->x-x3);
		y3-=4.*(scrk2->y-y3);
	}
    
    // second point, but if last point and it is exterior, extrapolate
	// forward 4 times the segment length
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
// The option at the end reverts to flat results in case want to proceed with that setting
// This method is no longer used, but can be interserted for debugging to verify hierachical crack crossing methods
void CrackHeader::CFFlatCrossing(double x1,double y1,double x2,double y2,Vector *norm,short *vfld,int pnum,int nodenum)
{
    Vector flatNorm;
    short svfld = *vfld;
    short flatVfld = FlatCrackCrossTest(x1,y1,x2,y2,&flatNorm);
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

// Determine if line from particle (x1,y1) to node (x2,y2) crosses this crack
// Return ABOVE_CRACK (1), BELOW_CRACK (2), or NO_CRACK (0) and crack normal in norm
// This method uses hierarchical crack in a binary tree
short CrackHeader::CrackCrossOnce(double x1,double y1,double x2,double y2,CrackSegment **crossSeg) const
{
    // recursive method to traverse tree hierarchy
    return CrackCrossLeafOnce(rootLeaf,x1,y1,x2,y2,crossSeg);
}

// Recursive Method to process each leaf in hierarchical traversal
short CrackHeader::CrackCrossLeafOnce(CrackLeaf *leaf,double x1,double y1,double x2,double y2,CrackSegment **crossSeg) const
{
    // check extents this leaf, return current cross if not in this leaf's extent
	if(!LineIsInExtents(x1,y1,x2,y2,leaf->cnear,leaf->cfar)) return NO_CRACK;
    
    // It is in extent of this leaf
    // if not terminal, go on to the child leaves
	short cross = NO_CRACK;
    if(!leaf->ChildrenAreSegments())
    {   CrackLeaf *child1,*child2;
        leaf->GetChildLeaves(&child1,&child2);
        cross = CrackCrossLeafOnce(child1,x1,y1,x2,y2,crossSeg);
		if(cross!=NO_CRACK) return cross;
        if(child2!=NULL) cross = CrackCrossLeafOnce(child2,x1,y1,x2,y2,crossSeg);
        return cross;
    }
    
    // Method 1: This code checks each segment now in a subroutine
	Vector norm;
    CrackSegment *scrk1,*scrk2;
    leaf->GetChildSegments(&scrk1,&scrk2);
    cross = CrackCrossOneSegment(scrk1,x1,y1,x2,y2,&norm,cross);
	if(cross!=NO_CRACK)
	{	*crossSeg = scrk1;
		return cross;
	}
    cross = CrackCrossOneSegment(scrk2,x1,y1,x2,y2,&norm,cross);
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
bool CrackHeader::CreateHierarchy(void)
{
    // all cracks have at least 2 points (and 1 segment)
    CrackSegment *scrk1 = firstSeg;
    if(scrk1==NULL) return false;
    CrackSegment *scrk2 = scrk1->nextSeg;
    if(scrk2==NULL) return false;
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
	if(!LineIsInExtents(x1,y1,x2,y2,scrk1->cnear,scrk1->cfar)) return NULL;

    // find intersection and see if in the contour
    double dx,dy;
    if(orient==HORIZONTAL)
    {   dy = scrk2->y-scrk1->y;
        if(DbleEqual(dy,0.)) return (x1==scrk1->x);		// parallel lines
        dx = scrk2->x-scrk1->x;
        crossPt->y = y1;
        crossPt->x = (crossPt->y-scrk1->y)*dx/dy + scrk1->x;
        if(x1<x2)
            return (crossPt->x>x1 && crossPt->x<=x2);
        else
            return (crossPt->x>x2 && crossPt->x<=x1);
    }
    
    // rest is VERTICAL
    dx = scrk2->x-scrk1->x;
    if(DbleEqual(dx,0.)) return (y1==scrk1->y);         // parallel lines
    dy = scrk2->y-scrk1->y;
    crossPt->x = x1;
    crossPt->y = (crossPt->x-scrk1->x)*dy/dx + scrk1->y;
    if(y1<y2)
        return (crossPt->y>y1 && crossPt->y<=y2);
    else
        return (crossPt->y>y2 && crossPt->y<=y1);
}

// If current J contour crossed any cracks check those cracks for crossing
// and if it crosses, stop the grow at that crack
// Return relative change made in grow
double CrackHeader::AdjustGrowForCrossing(Vector *grow,CrackSegment *crkTip)
{
	// exit if no in the contour
	if(crossedCracks == NULL) return 1.0;
	
	// get growing line segment
	double x1 = crkTip->x;
	double y1 = crkTip->y;
	double x2 = x1 + grow->x;
	double y2 = y1 + grow->y;
	
	// loop over crosse cracks
	double p = 1.0;
	CrackSegment *crossSeg;
	CrossedCrack *nextCross = (CrossedCrack *)crossedCracks->firstObject;
	while(nextCross!=NULL)
	{	CrackHeader *nextCrack = nextCross->crack;
		if(nextCrack->CrackCrossOnce(x1,y1,x2,y2,&crossSeg)!=NO_CRACK)
		{	// propagation path is x = x1 + p*grow->x, y = y1 + p*grow->y
			// line along the crack segment is x = crossSeg->x + s*dxs, y = crossSeg->y + s*dys
			CrackSegment *nextSeg = crossSeg->nextSeg;
			double dxs = nextSeg->x-crossSeg->x;
			double dys = nextSeg->y-crossSeg->y;
			
			// equate x and solve for s : s = (x1 - crossSeg->x + p*grow->x)/dxs
			//    But, if dxs==0 then x1 + p*grow->x = crossSeg->x to give p = (crossSeg->x-x1)/grow->x
			//		   if grow->x==0 too, parallel lines so do nothing
			// equate y and solve for p : p*(grow->y*dxs - grow->x*dys) = (crossSeg->y - y1)*dxs + (x1 - crossSeg->x)*dys
			//    But, if (grow->y*dxs - grow->x*dys)==0, parallel lines so do nothing
			double cp = grow->y*dxs - grow->x*dys;			// will be zero if lines are parallel
			if(cp!=0.)
			{	if(dxs==0.)
					p = (crossSeg->x-x1)/grow->x;
				else
					p = ((crossSeg->y - y1)*dxs + (x1 - crossSeg->x)*dys)/cp;
				
				// check last segment length
				CrackSegment *prevSeg = crkTip->nextSeg;
				if(prevSeg==NULL) prevSeg = crkTip->prevSeg;
				double pxs = crkTip->x-prevSeg->x;
				double pys = crkTip->y-prevSeg->y;
				double plen = sqrt(pxs*pxs+pys*pys);					// length of crack tip segment
				double glen = sqrt(grow->x*grow->x+grow->y*grow->y);	// length to proposed new tip segment
				if(PropagateTask::cellsPerPropagationStep>.7)
				{	int numSegs= 2*(PropagateTask::cellsPerPropagationStep+.25);
					glen /= (double)numSegs;
				}
				double pmin = plen/glen<0.5 ? 0.5 : 0.0;
				
				// Adjust grow, but do not put put two small segments in a row
				if(p>=pmin && p<=1.0)
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

#pragma mark ACCESSORS

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
int CrackHeader::CriterionNeeds(bool &totalEnergyNeeded)
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
		thisCrackNeeds|=theMaterials[tipCrk->tipMatnum-1]->CriterionNeeds(0,totalEnergyNeeded);
		if(GetAllowAlternate(crkTipIdx))
			thisCrackNeeds|=theMaterials[tipCrk->tipMatnum-1]->CriterionNeeds(1,totalEnergyNeeded);
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

// Check is lines is within extent of leaf or segment
// return false if line within extents or false otherwise
// See JANOSU-6-66
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
