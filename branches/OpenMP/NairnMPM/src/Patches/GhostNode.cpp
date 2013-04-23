/********************************************************************************
	GridNode.cpp
	NairnMPM

	Created by John Nairn on 4/11/13.
	Copyright (c) 2013 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Patches/GhostNode.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Nodes/NodalPoint2D.hpp"
#include "Nodes/NodalPoint3D.hpp"
#include "Custom_Tasks/TransportTask.hpp"

#pragma mark INITIALIZATION

// Constructors
// Create for node at zero-based (row,col) grid position
GhostNode::GhostNode(int row,int col,bool interiorRow,bool interiorCol)
{
	// create ghost node for the real node
	int xbase = row*mpmgrid.yplane;			// node at start of this row
	int realNode = xbase+col+1;
	
	// get the real node
	if(row>=0 && col>=0 && !mpmgrid.EdgeNode(row,'y') && !mpmgrid.EdgeNode(col,'x'))
		real = nd[realNode];
	else
		real = NULL;
	
	// should we create a ghost node
	bool makeGhost = real!=NULL;
	if(makeGhost)
	{	if(mpmgrid.EdgeNode(row+1,'y'))
		{	// top row - no ghost if top-right of grid or interior column
			if(interiorCol || mpmgrid.EdgeNode(col+1,'x'))
				makeGhost = FALSE;
		}
		else if(mpmgrid.EdgeNode(col+1,'x') && interiorRow)
		{	// right edge or interior row needs no ghost
			makeGhost = FALSE;
		}
	}
	
	// create ghost node now
	if(makeGhost)
		ghost = new NodalPoint2D(real);
	else
		ghost = NULL;
		
	//cout << "...... node = " << realNode << "," << ghost << "," << real << "," << (real==ghost) << endl;
}

// Create for node at zero-based (row,col) grid position
GhostNode::GhostNode(int row,int col,int rank,bool interiorRow,bool interiorCol,bool interiorRank)
{
	// create ghost node for the real node
	int xbase = rank*mpmgrid.zplane + row*mpmgrid.yplane;		// node at start of this row
	int realNode = xbase+col+1;
	
	// get the real node
	if(row>=0 && col>=0 && rank>=0 && !mpmgrid.EdgeNode(row,'y') && !mpmgrid.EdgeNode(col,'x') && !mpmgrid.EdgeNode(rank,'z'))
		real = nd[realNode];
	else
		real = NULL;
	
	// should we create a ghost node
	bool makeGhost = real!=NULL;
	if(makeGhost)
	{	if(mpmgrid.EdgeNode(rank+1,'z'))
		{	// apex rank
			if(interiorRow && interiorCol)
			{	// no ghost for interior node on rank edge
				makeGhost = FALSE;
			}
			else if(mpmgrid.EdgeNode(row+1,'y'))
			{	// no ghost on y edge or on top-right edge of grid
				if(interiorCol || mpmgrid.EdgeNode(col+1,'x'))
					makeGhost = FALSE;
			}
			else if(mpmgrid.EdgeNode(col+1,'x') && interiorRow)
			{	// no ghost on x edge (top-right handled above)
				makeGhost = FALSE;
			}
		}
		else if(interiorRank)
		{	if(mpmgrid.EdgeNode(row+1,'y'))
			{	// top row - no ghost if interior column
				if(interiorCol || mpmgrid.EdgeNode(col+1,'x'))
					makeGhost = FALSE;
			}
			else if(mpmgrid.EdgeNode(col+1,'x') && interiorRow)
			{	// right edge or interior row needs no ghost
				makeGhost = FALSE;
			}
		}
	}
	
	// create ghost node now
	if(makeGhost)
		ghost = new NodalPoint3D(real);
	else
		ghost = NULL;
	
	//cout << "...... node = " << realNode << "," << ghost << "," << real << "," << (real==ghost) << endl;
}


#pragma mark GhostNode: Methods

// initialize ghost nodes for next time step
void GhostNode::InitializeForTimeStep()
{	if(ghost!=NULL)
		ghost->InitializeForTimeStep();
}

// copy allocated velocity fields from real node to ghost node
void GhostNode::InitializationReduction(void)
{	if(ghost!=NULL)
		real->CopyFieldInitialization(ghost);
}

// When Grid Forces task is done transfer ghost node force to real nodes
void GhostNode::MassAndMomentumReduction(void)
{	if(ghost!=NULL)
	{	ghost->CopyMassAndMomentum(real);
		
		// transport calculations
		TransportTask *nextTransport=transportTasks;
		while(nextTransport!=NULL)
			nextTransport=nextTransport->Task1Reduction(real,ghost);
	}
}

// initialize ghost nodes for next time step
void GhostNode::RezeroNodeTask6(double delTime)
{	if(ghost!=NULL)
        ghost->RezeroNodeTask6(delTime);
}

// When Grid Forces task is done transfer ghost node force to real nodes
void GhostNode::MassAndMomentumReductionLast(void)
{	if(ghost!=NULL)
        ghost->CopyMassAndMomentumLast(real);
}

// When Grid Forces task is done transfer ghost node force to real nodes
void GhostNode::GridForcesReduction(void)
{	if(ghost!=NULL)
	{	ghost->CopyGridForces(real);
		
		// transport forces
		TransportTask *nextTransport=transportTasks;
		while(nextTransport!=NULL)
			nextTransport=nextTransport->CopyForces(real,ghost);
	}
}

// initialize ghost nodes for next time step
void GhostNode::ZeroDisp(void)
{	if(ghost!=NULL)
        ghost->ZeroDisp(real);
}

// When Grid Forces task is done transfer ghost node force to real nodes
void GhostNode::JKTaskReduction(void)
{	if(ghost!=NULL)
        ghost->CopyUGradientStressEnergy(real);
}

// initialize ghost nodes for next time step
void GhostNode::DeleteDisp(void)
{	if(ghost!=NULL)
        ghost->DeleteDisp(real);
}

#pragma mark GhostNode: Accessors

// if has ghost node return it, otherwise return real node
NodalPoint *GhostNode::GetNodePointer(void)
{	if(ghost==NULL && real==NULL)
		throw "double NULL ghost node";
	return ghost!=NULL ? ghost : real ;
}

