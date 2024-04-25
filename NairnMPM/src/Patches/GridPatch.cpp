/********************************************************************************
    GridPatch.cpp
    nairn-mpm-fea

    Created by John Nairn on 4/10/13.
    Copyright (c) 2013 John A. Nairn, All rights reserved.
 
	The ghost nodes wrap around the element
	In 2D (on interior slices in 3D) have
		first ghostRows rows: xnodes + 2*ghostRows
		next ynodes-1 rows: 2*ghostRows+1
		last ghostRows+1 rows: xnodes + 2*ghostRows
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Patches/GridPatch.hpp"
#include "Patches/GhostNode.hpp"
#include "Nodes/NodalPoint.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "NairnMPM_Class/Reservoir.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "Exceptions/CommonException.hpp"

// globals
GridPatch **patches;            // list of patches (or NULL if only one patch or if serial)
int GridPatch::ghostRows = 1;   // number of ghost rows. If needed, increase for higher strain limits

#pragma mark GridPath: Initialization

// Constructors
GridPatch::GridPatch(int xmin,int xmax,int ymin,int ymax,int zmin,int zmax)
{	
    x0 = xmin-1;
    x1 = xmax-1;
    y0 = ymin-1;
    y1 = ymax-1;
    z0 = zmin-1;			  // = 0 if 2D patch
    z1 = zmax-1;              // = -1 if 2D patch
    //cout << " (" << x0+1 << "-" << x1+1 << "),(" << y0+1 << "-" << y1+1 << "),("  << z0+1 << "-" << z1+1 << ")" << endl;
	//cout << "----------------------------------------------------------------" << endl;
	
	xn = x1-x0+1;		// number of elements in x direction (xn+1 nodes)
	yn = y1-y0+1;		// number of elements in y direction (yn+1 nodes)
	zn = z1-z0+1;		// number of elements in z direction (zn+1 nodes) (=0 if 2D)
	
	// particle pointers
	firstNR = NULL;
	firstRB = NULL;
	firstRC = NULL;
	firstRBC = NULL;
	
	ghosts = NULL;
	numGhosts=0;
	
	lastToMove = NULL;
}

// Create all ghost nodes
// throws std::bad_alloc
bool GridPatch::CreateGhostNodes(void)
{
	int i,j,k;
	
	// 2D Patches:
	//   Each patch owns at x0 and y0, but nodes at x1 and y1 belong to next patch (or are edge of the grid)
	//   Need ghostRows on left and right and one at x1 for interior rows with interiorRow=(2*ghostRows+1) nodes
	//   A full row of nodes counting ghost nodes as fullRow=(xn+1)+2*ghostRows=xn+interiorRow nodes
	//   basePartial is zero based address of first interior row
	//   baseTop is zero based address of first full row on the top (have ghostRows+1 of these rows
	//   Total ghost nodes = (2*ghostRows+1)*fullRow+yn*interiorRow = interiorRow*(fullRow+yn)
	
	// partial counts
	interiorRow = 2*ghostRows+1;					// midrow count (2D patch) AND # of bottom/base+top/apex rows
	fullRow = xn + interiorRow;						// top and bottom row counts, 2D patch
	basePartial = ghostRows*fullRow;				// start of mid rows, 2D patch (0 based)
	baseTop = basePartial + yn*interiorRow;			// start of top full rows, 2D patch (0 based)
    //cout << "... base partial " << basePartial << ", base top " << baseTop << endl;
	
	// create ghost nodes
	if(zn<=0)
	{	// 2D = one rank
		numGhosts = interiorRow*(fullRow+yn);
	}
	else
	{	// 3D ranks
        interiorRank = interiorRow*(fullRow+yn);		// midblock rank count, 3D patch (same as total in 2D patch)
        fullRank = fullRow*(yn + interiorRow);			// base and apex rank counts, 3D patch
        baseInterior = ghostRows*fullRank;              // start of interior ranks, 3D patch (0 based)
        baseApex = baseInterior + zn*interiorRank;      // start of top full ranks, 3D patch (0 based)
		numGhosts = interiorRow*fullRank + zn*interiorRank;
        //cout << "... base interior " << baseInterior << ", base apex " << baseApex << endl;
	}
	ghosts = new (nothrow) GhostNode *[numGhosts];
	if(ghosts==NULL) return false;
	//cout << "... gtot = " << numGhosts << endl;

	if(zn<=0)
	{	// 2D row by row
		int row,col,g=0;
		for(i=-ghostRows;i<=yn+ghostRows;i++)
		{	row = y0+i;					// 0 based grid row
			//cout << endl;
			for(j=-ghostRows;j<=xn+ghostRows;j++)
			{	// skip interior nodes
				if(i>=0 && i<yn && j>=0 && j<xn) continue;
				col = x0+j;				// 0 based gid row
				//cout << "...(" << g << "," << row << "," << col << ")" << endl;
				ghosts[g] = new GhostNode(row,col,i>=0 && i<yn,j>=0 && j<xn);
				if(ghosts[g]==NULL) return false;
				g++;
			}
		}
	}
	else
	{	// rank by row by row
		int row,col,rank,g=0;
		for(k=-ghostRows;k<=zn+ghostRows;k++)
		{	rank = z0+k;					// 0 based grid rank
			//cout << "\n**** begin rank on z rank = " << rank << endl;
			for(i=-ghostRows;i<=yn+ghostRows;i++)
			{	row = y0+i;					// 0 based grid row
				//cout << endl;
				for(j=-ghostRows;j<=xn+ghostRows;j++)
				{	// skip interior nodes
					if(i>=0 && i<yn && j>=0 && j<xn && k>=0 && k<zn) continue;
					col = x0+j;				// 0 based gid row
					//cout << "...(" << g << "," << row << "," << col << "," << rank << ")" << endl;
					ghosts[g] = new GhostNode(row,col,rank,i>=0 && i<yn,j>=0 && j<xn,k>=0 && k<zn);
					if(ghosts[g]==NULL) return false;
					g++;
				}
			}
		}
	}

	return true;
}

#pragma mark GridPatch: Methods

// initialize ghost nodes for next time step
void GridPatch::InitializeForTimeStep(void)
{   for(int i=0;i<numGhosts;i++)
		ghosts[i]->InitializeForTimeStep();
}

// When initialization is done, copy velocity field to real nodes
void GridPatch::InitializationReduction(void)
{	for(int i=0;i<numGhosts;i++)
		ghosts[i]->InitializationReduction();
}

// When Mass and momentum task is done transfer ghost node values to real nodes
void GridPatch::MassAndMomentumReduction(void)
{	for(int i=0;i<numGhosts;i++)
		ghosts[i]->MassAndMomentumReduction();
}

// Support XPIC calculations (always xpicCalculation==COPY_VSTARNEXT, and k-pass in XPIC loop)
void GridPatch::XPICSupport(int xpicCalculation,int xpicOption,NodalPoint *real,double timestep,int m,int k,double vsign)
{	for(int i=0;i<numGhosts;i++)
	ghosts[i]->XPICSupport(xpicCalculation,xpicOption,real,timestep,m,k,vsign);
}

// zero gTnext on ghost nodes
void GridPatch::InitializeForXPICTransport(TransportTask *nextTransport,int xpicOption)
{	for(int i=0;i<numGhosts;i++)
		ghosts[i]->InitializeForXPICTransport(nextTransport,xpicOption);
}

void GridPatch::XPICReductionTransport(int k,TransportTask *nextTransport)
{	for(int i=0;i<numGhosts;i++)
		ghosts[i]->XPICReductionTransport(k,nextTransport);
}

// initialize ghost nodes for next time step
void GridPatch::RezeroNodeTask6(double delTime)
{   for(int i=0;i<numGhosts;i++)
        ghosts[i]->RezeroNodeTask6(delTime);
}

// When Grid Forces task is done transfer ghost node force to real nodes
void GridPatch::MassAndMomentumReductionLast(void)
{	for(int i=0;i<numGhosts;i++)
        ghosts[i]->MassAndMomentumReductionLast();
}

// When Grid Forces task is done transfer ghost node force to real nodes
void GridPatch::GridForcesReduction(void)
{	for(int i=0;i<numGhosts;i++)
		ghosts[i]->GridForcesReduction();
}

// initialize ghost nodes for next time step
void GridPatch::ZeroDisp()
{   for(int i=0;i<numGhosts;i++)
        ghosts[i]->ZeroDisp();
}

// When Grid Forces task is done transfer ghost node force to real nodes
void GridPatch::JKTaskReduction(void)
{	for(int i=0;i<numGhosts;i++)
        ghosts[i]->JKTaskReduction();
}

// initialize ghost nodes for next time step
void GridPatch::DeleteDisp(void)
{   for(int i=0;i<numGhosts;i++)
        ghosts[i]->DeleteDisp();
}

// prepare particle to be moved to another patch
bool GridPatch::AddMovingParticle(MPMBase *mptr,GridPatch *newPatch,MPMBase *prevMptr)
{	// create data structure
	MovingData *nextToMove = new MovingData;
	
	// fill with data
	nextToMove->movingMptr = mptr;
	nextToMove->newPatch = (void *)newPatch;
	nextToMove->previousMptr = prevMptr;
	nextToMove->previousMoveData = (void *)lastToMove;
	
	// reset last to move and exit
	lastToMove = nextToMove;
	return true;
}

// If any scheduled to move, move them now and clear data in the process
void GridPatch::MoveParticlesToNewPatches(void)
{	// move if any need to be moved
	while(lastToMove!=NULL)
	{	// move to new one or delete if new one is NULL
		GridPatch *patch = (GridPatch *)lastToMove->newPatch;
		if(patch!=NULL)
		{	// remove from previous patch
			RemoveParticleAfter(lastToMove->movingMptr,lastToMove->previousMptr);
			// add to new one
			patch->AddParticle(lastToMove->movingMptr);
		}
		else
			mpmReservoir->DeleteParticle(lastToMove->movingMptr);
		
		// find next one and delete old one
		MovingData *nextToMove = (MovingData *)lastToMove->previousMoveData;
		delete lastToMove;
		
		// on to next one
		lastToMove = nextToMove;
	}
}

// add particle to this patch by
//    1. Set is next particle the previous first of that type
//    2. Set first of this type to the added particle
void GridPatch::AddParticle(MPMBase *mptr)
{
	const MaterialBase *matref = theMaterials[mptr->MatID()];		// material object for this particle
	
	if(matref->IsRigid())
	{	// sf it BC, contact, or block
		if(matref->IsRigidBC())
		{	mptr->SetNextObject(firstRBC);
			firstRBC = mptr;
		}
		else if(matref->IsRigidContact())
		{	mptr->SetNextObject(firstRC);
			firstRC = mptr;
		}
		else
		{	mptr->SetNextObject(firstRB);
			firstRB = mptr;
		}
	}
	else
	{	mptr->SetNextObject(firstNR);
		firstNR = mptr;
	}
}

// Remove the particle mptr by finding prevPtr and then use standard remove code
// If you know the previous particle, calling RemoveParticleAfter() is more efficient
void GridPatch::RemoveParticle(MPMBase *mpmptr)
{
	for(int block=FIRST_NONRIGID;block<=FIRST_RIGID_BC;block++)
	{	// get first material point in this block
		MPMBase *mptr = GetFirstBlockPointer(block);
		MPMBase *prevMptr = NULL;
		while(mptr!=NULL)
		{	if(mptr==mpmptr)
			{	RemoveParticleAfter(mpmptr,prevMptr);
				return;
			}
			prevMptr = mptr;
			mptr = (MPMBase *)prevMptr->GetNextObject();
		}
	}
}

// Remove the particle mptr from this patch that is after prevMptr
void GridPatch::RemoveParticleAfter(MPMBase *mptr,MPMBase *prevMptr)
{
	// material point after the one being removed (might be NULL)
	MPMBase *nextMptr = (MPMBase *)mptr->GetNextObject();
	
	// If prevPtr is NULL, then the one being removed is the first of
	// a certain type. Need to reset the first pointer
	if(prevMptr==NULL)
	{	// Check material type of particle being removed
		const MaterialBase *matref = theMaterials[mptr->MatID()];		// material object for this particle
		if(matref->IsRigid())
		{	if(matref->IsRigidBC())
			{	if(mptr == firstRBC)
				{	// still first, switch to next one as new first
					firstRBC = nextMptr;
					return;
				}
				prevMptr = firstRBC;
			}
			else if(matref->IsRigid())
			{	if(mptr == firstRC)
				{	// still first, switch to next one as new first
					firstRC = nextMptr;
					return;
				}
				prevMptr = firstRC;
			}
			else
			{	if(mptr == firstRB)
				{	// still first, switch to next one as new first
					firstRB = nextMptr;
					return;
				}
				prevMptr = firstRB;
			}
		}
		else
		{	if(mptr == firstNR)
			{	// still first, switch to next one as new first
				firstNR = nextMptr;
				return;
			}
			prevMptr = firstNR;
		}
		
		// I think above code always exits, if not, this finds pointer in appropriate list
		// prevPtr is now first of one class
		// step to previous first particle to find its new previous point
		MPMBase *currentMptr = (MPMBase *)prevMptr->GetNextObject();
		while(currentMptr != mptr)
		{	prevMptr = currentMptr;
			currentMptr = (MPMBase *)prevMptr->GetNextObject();
		}
	}
	
	// removing one in the middle or at the end of linked list
	prevMptr->SetNextObject(nextMptr);

}

#pragma mark GridPatch: Accessors

// return pointer to real or ghost node for 1-based node number num in the global grid
//	as appropriate for this patch
// throws CommonException()
NodalPoint *GridPatch::GetNodePointer(int num)
{
	// if single patch, use the real node
	if(ghosts==NULL) return nd[num];
	
	// look for ghost node
	int g,col,row;
	
	// 2D
	if(zn<=0)
	{	// get rwo column in this patch
		col = (num-1)%mpmgrid.yplane;               // row, col in global grid
		row = (num-1)/mpmgrid.yplane;
		col -= x0;									// zero based within the patch
		row -= y0;
		
		// is it an owned node?
		if(row>=0 && row<yn && col>=0 && col<xn) return nd[num];
		
		// is it out of this patch
		if(row<-ghostRows || row>yn+ghostRows || col<-ghostRows || col>xn+ghostRows)
			throw CommonException("Need ghost node that is outside this patch (i.e., increase ghost rows)","GridPatch::GetNodePointer");
		
		if(row<0)
		{	// ghost in full rows near the bottom
			g = (row+ghostRows)*fullRow + ghostRows + col;
		}
		
		else if(row<yn)
		{	// ghost within partial rows
			if(col<0)
				g = basePartial + row*interiorRow + ghostRows + col ;
			else
				g = basePartial + row*interiorRow + ghostRows + col - xn;
		}
		
		else
		{	// top full rows
			g = baseTop + (row-yn)*fullRow + ghostRows + col;
		}
	}
	
	else
	{	// 3D patch
        int rank = (num-1)/mpmgrid.zplane;          // 0-based rank in global grid
        int rnum = (num-1)%mpmgrid.zplane;          // 0-based number in rank
		col = rnum%mpmgrid.yplane;					// 0-baed row, col in global grid
		row = rnum/mpmgrid.yplane;
		col -= x0;									// zero based within the patch
		row -= y0;
        rank -= z0;
        
		// is it an owned node?
		if(row>=0 && row<yn && col>=0 && col<xn && rank>=0 && rank<zn) return nd[num];
        
		// is it out of this patch
		if(row<-ghostRows || row>yn+ghostRows || col<-ghostRows || col>xn+ghostRows || rank<-ghostRows || rank>zn+ghostRows)
        {   cout << "# ghost for node " << num << " to (" << row << "," << col << "," << rank << ") outside patch" << endl;
			throw CommonException("Need ghost node that is outside this patch (i.e., increase ghost rows)","GridPatch::GetNodePointer");
        }
        
        if(rank<0)
        {   // ghost in full slices below base of the element
            g = (rank+ghostRows)*fullRank + (row+ghostRows)*fullRow + ghostRows + col;
        }
        else if(rank<zn)
        {   // ghost within partial slices
            int baseRank = baseInterior + rank*interiorRank;
            if(row<0)
            {	// ghost in full rows near the bottom
                g = baseRank + (row+ghostRows)*fullRow + ghostRows + col;
            }
            
            else if(row<yn)
            {	// ghost within partial rows
                if(col<0)
                    g = baseRank + basePartial + row*interiorRow + ghostRows + col ;
                else
                    g = baseRank + basePartial + row*interiorRow + ghostRows + col - xn;
            }
            
            else
            {	// top full rows
                g = baseRank + baseTop + (row-yn)*fullRow + ghostRows + col;
            }
        }
        
        else
        {   // ghost in full slices above the patch
            g = baseApex + (rank-zn)*fullRank + (row+ghostRows)*fullRow + ghostRows + col;
        }
	}

	// if ghosts[g] has ghost node return it, otherwise return real node
	// a ghosts[g] with no nodes should not reach here
	if(g<0 || g>=numGhosts)
	{	cout << "# ghost for node " << num << " out of range (" << g << ") for (" << row << "," << col << ")"
				<< ") from (" << y0 << "," << x0 << ")" << endl;
		throw CommonException("ghost index out of range","GridPatch::GetNodePointer");
	}
	NodalPoint *thePtr = ghosts[g]->GetNodePointer();
	if(thePtr==NULL)
	{	cout << "NULL pointer for " << g << " node " << num << " at (" << row << "," << col
				<< ") from (" << y0 << "," << x0 << ")" << endl;
		throw CommonException("NULL node pointer","GridPatch::GetNodePointer");
	}
	
	return thePtr;
}

// return pointer to real or ghost node for 1=based node number num in the global grid
//	as appropriate for this patch
// throws CommonException()
NodalPoint *GridPatch::GetNodePointer(int num,bool debug)
{
	// if single patch, use the real node
	if(ghosts==NULL) return nd[num];
	
	// look for ghost node
	int g,col,row;
	
	// 2D
	if(zn<=0)
	{	// get rwo column in this patch
		col = (num-1)%mpmgrid.yplane;               // row, col in global grid
		row = (num-1)/mpmgrid.yplane;
		col -= x0;									// zero based within the patch
		row -= y0;
		
		// is it an owned node?
		if(row>=0 && row<yn && col>=0 && col<xn) return nd[num];
		
		// is it out of this patch
		if(row<-ghostRows || row>yn+ghostRows || col<-ghostRows || col>xn+ghostRows)
			throw CommonException("Need ghost node that is outside this patch (i.e., increase ghost rows)","GridPatch::GetNodePointer");
		
		if(row<0)
		{	// ghost in full rows near the bottom
			g = (row+ghostRows)*fullRow + ghostRows + col;
		}
		
		else if(row<yn)
		{	// ghost within partial rows
			if(col<0)
				g = basePartial + row*interiorRow + ghostRows + col ;
			else
				g = basePartial + row*interiorRow + ghostRows + col - xn;
		}
		
		else
		{	// top full rows
			g = baseTop + (row-yn)*fullRow + ghostRows + col;
		}
	}
	
	else
	{	// 3D patch
        int rank = (num-1)/mpmgrid.zplane;          // 0-based rank in global grid
        int rnum = (num-1)%mpmgrid.zplane;          // 0-based number in rank
		col = rnum%mpmgrid.yplane;					// 0-baed row, col in global grid
		row = rnum/mpmgrid.yplane;
        cout << "row,col,rank = " << row << "," << col << "," << rank << endl;
		col -= x0;									// zero based within the patch
		row -= y0;
        rank -= z0;
        cout << "patch row,col,rank = " << row << "," << col << "," << rank << endl;
        
		// is it an owned node?
		if(row>=0 && row<yn && col>=0 && col<xn && rank>=0 && rank<zn) return nd[num];
        
		// is it out of this patch
		if(row<-ghostRows || row>yn+ghostRows || col<-ghostRows || col>xn+ghostRows || rank<-ghostRows || rank>zn+ghostRows)
        {   cout << "# ghost for node " << num << " to (" << row << "," << col << "," << rank << ") outside patch" << endl;
			throw CommonException("Need ghost node that is outside this patch (i.e., increase ghost rows)","GridPatch::GetNodePointer");
        }
        
        if(rank<0)
        {   // ghost in full slices below base of the element
            g = (rank+ghostRows)*fullRank + (row+ghostRows)*fullRow + ghostRows + col;
        }
        else if(rank<zn)
        {   // ghost within partial slices
            int baseRank = baseInterior + rank*interiorRank;
            if(row<0)
            {	// ghost in full rows near the bottom
                g = baseRank + (row+ghostRows)*fullRow + ghostRows + col;
            }
            
            else if(row<yn)
            {	// ghost within partial rows
                if(col<0)
                    g = baseRank + basePartial + row*interiorRow + ghostRows + col ;
                else
                    g = baseRank + basePartial + row*interiorRow + ghostRows + col - xn;
            }
            
            else
            {	// top full rows
                g = baseRank + baseTop + (row-yn)*fullRow + ghostRows + col;
            }
        }
        
        else
        {   // ghost in full slices above the patch
            g = baseApex + (rank-zn)*fullRank + (row+ghostRows)*fullRow + ghostRows + col;
        }
	}
    
	// if ghosts[g] has ghost node return it, otherwise return real node
	// a ghosts[g] with no nodes should not reach here
	if(g<0 || g>=numGhosts)
	{	cout << "ghost for node " << num << " out of range (" << g << ") for (" << row << "," << col << ")"
        << ") from (" << y0 << "," << x0 << ")" << endl;
		throw CommonException("ghost index out of range","GridPatch::GetNodePointer");
	}
	NodalPoint *thePtr = ghosts[g]->GetNodePointer();
	if(thePtr==NULL)
	{	cout << "NULL pointer for " << g << " node " << num << " at (" << row << "," << col
        << ") from (" << y0 << "," << x0 << ")" << endl;
		throw CommonException("NULL node pointer","GridPatch::GetNodePointer");
	}
	
	return thePtr;
}


// return material point pointer by ID 0 to 3 for FIRST_NONRIGID=0,
// FIRST_RIGID_BLOCK, FIRST_RIGID_CONTACT, FIRST_RIGID_BC
MPMBase *GridPatch::GetFirstBlockPointer(int block)
{
	switch(block)
	{	case FIRST_NONRIGID:
			return firstNR;
		case FIRST_RIGID_CONTACT:
			return firstRC;
		case FIRST_RIGID_BC:
			return firstRBC;
		case FIRST_RIGID_BLOCK:
			return firstRB;
		default:
			break;
	}
	return NULL;
}

// return private variables
GhostNode **GridPatch::GetGhosts(int *getNum)
{	*getNum = numGhosts;
	return ghosts;
}
