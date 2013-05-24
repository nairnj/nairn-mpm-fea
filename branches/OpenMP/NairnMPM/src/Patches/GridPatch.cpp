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

#include "Patches/GridPatch.hpp"
#include "Patches/GhostNode.hpp"
#include "Nodes/NodalPoint.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/MaterialBase.hpp"

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
	zn = z1-z0+1;		// number of elements in z direction (zn+1 nodes) (=0 if 3D)
	
	// particle pointers
	firstNR = NULL;
	firstRC = NULL;
	firstRBC = NULL;
	
	ghosts = NULL;
	numGhosts=0;
}

// Create all ghost nodes
bool GridPatch::CreateGhostNodes(void)
{
	int i,j,k;
	
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
	ghosts = (GhostNode **)malloc(numGhosts*sizeof(GhostNode));
	if(ghosts==NULL) return FALSE;
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
				if(ghosts[g]==NULL) return FALSE;
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
					if(ghosts[g]==NULL) return FALSE;
					g++;
				}
			}
		}
	}

	return TRUE;
}

#pragma mark GridPatch: Methods

// initialize ghost nodes for next time step
void GridPatch::InitializeForTimeStep()
{   for(int i=0;i<numGhosts;i++)
		ghosts[i]->InitializeForTimeStep();
}

// When initialization is done, copy velocity field to real nodes
void GridPatch::InitializationReduction(void)
{	for(int i=0;i<numGhosts;i++)
		ghosts[i]->InitializationReduction();
}

// When Grid Forces task is done transfer ghost node force to real nodes
void GridPatch::MassAndMomentumReduction(void)
{	for(int i=0;i<numGhosts;i++)
		ghosts[i]->MassAndMomentumReduction();
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

// add particle to this patch - it is put at the beginning
void GridPatch::AddParticle(MPMBase *mptr)
{
	const MaterialBase *matref = theMaterials[mptr->MatID()];		// material object for this particle
	
	if(matref->RigidBC())
	{	mptr->SetNextObject(firstRBC);
		firstRBC = mptr;
	}
	else if(matref->Rigid())
	{	mptr->SetNextObject(firstRC);
		firstRC = mptr;
	}
	else
	{	mptr->SetNextObject(firstNR);
		firstNR = mptr;
	}
}

// add particle to this patch - it is put at the beginning
// in block 0 to 2 (NR, RC, RBC)
void GridPatch::RemoveParticleAfter(MPMBase *mptr,MPMBase *prevMptr)
{
	// material point after the one being removed
	MPMBase *nextMptr = (MPMBase *)mptr->GetNextObject();
	
	if(prevMptr!=NULL)
	{	// removing one in the middle or at the end of linked list
		prevMptr->SetNextObject(nextMptr);
	}
	else
	{	// removing the first particle
		const MaterialBase *matref = theMaterials[mptr->MatID()];		// material object for this particle
		if(matref->RigidBC())
			firstRBC = nextMptr;
		else if(matref->Rigid())
			firstRC = nextMptr;
		else
			firstNR = nextMptr;
	}
}

#pragma mark GridPatch: Accessors

// return pointer to real or ghost node for 1=based node number num in the global grid
//	as appropriate for this patch
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
			throw "Need ghost node that is outside this patch";
		
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
        {   cout << num << " to (" << row << "," << col << "," << rank << ")" << endl;
			throw "Need ghost node that is outside this patch";
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
		throw "ghost index out of range";
	}
	NodalPoint *thePtr = ghosts[g]->GetNodePointer();
	if(thePtr==NULL)
	{	cout << "NULL pointer for " << g << " node " << num << " at (" << row << "," << col
				<< ") from (" << y0 << "," << x0 << ")" << endl;
		throw "NULL node pointer";
	}
	
	return thePtr;
}

// return pointer to real or ghost node for 1=based node number num in the global grid
//	as appropriate for this patch
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
			throw "Need ghost node that is outside this patch";
		
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
        {   cout << num << " to (" << row << "," << col << "," << rank << ")" << endl;
			throw "Need ghost node that is outside this patch";
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
		throw "ghost index out of range";
	}
	NodalPoint *thePtr = ghosts[g]->GetNodePointer();
	if(thePtr==NULL)
	{	cout << "NULL pointer for " << g << " node " << num << " at (" << row << "," << col
        << ") from (" << y0 << "," << x0 << ")" << endl;
		throw "NULL node pointer";
	}
	
	return thePtr;
}


// return material point pointer by ID 0 to 2  or FIRST_NONRIGID=0,FIRST_RIGID_CONTACT,FIRST_RIGID_BC
MPMBase *GridPatch::GetFirstBlockPointer(int block)
{
	switch(block)
	{	case FIRST_NONRIGID:
			return firstNR;
		case FIRST_RIGID_CONTACT:
			return firstRC;
		case FIRST_RIGID_BC:
			return firstRBC;
		default:
			break;
	}
	return NULL;
}

