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

#pragma mark INITIALIZATION

// Constructors
// Create for node at zero-based (row,col) grid position
GhostNode::GhostNode(int row,int col,bool interiorRow,bool interiorCol,bool is2D)
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
		{	// top row - no if top-right of grid or interior column
			if(mpmgrid.EdgeNode(col+1,'x') || interiorCol)
				makeGhost = FALSE;
		}
		else if(mpmgrid.EdgeNode(col+1,'x') && interiorRow)
		{	// right edge - no unless top of ths patch
			makeGhost = FALSE;
		}
	}
	
	// create ghost node now
	if(makeGhost)
	{	if(is2D)
			ghost = new NodalPoint2D(real);
		else
			ghost = new NodalPoint3D(real);
	}
	else
		ghost = NULL;
		
	cout << "...... node = " << realNode << "," << ghost << "," << real << "," << (real==ghost) << endl;
}

#pragma mark GhostNode: Methods

// initialize real and ghost nodes for next time step
void GhostNode::InitializeForTimeStep()
{
	if(ghost!=NULL)
	{	// If has ghost, initialize it (real is done by owner)
		ghost->InitializeForTimeStep();
	}
	else if(real!=NULL)
	{	// if no ghost, but a real it is edge node and needs initialization
		real->InitializeForTimeStep();
	}
}

// copy allocated velocity fields to real nodes
void GhostNode::InitializationReduction(void)
{	if(ghost!=NULL)
		ghost->CopyFieldInitialization(real);
}

// When Grid Forces task is done transfer ghost node force to real nodes
void GhostNode::MassAndMomentumReduction(void)
{	if(ghost!=NULL)
		ghost->CopyMassAndMomentum(real);
}

// When Grid Forces task is done transfer ghost node force to real nodes
void GhostNode::GridForcesReduction(void)
{	if(ghost!=NULL)
		ghost->CopyGridForces(real);
}

#pragma mark GhostNode: Accessors

// if has ghost node return it, otherwise return real node
NodalPoint *GhostNode::GetNodePointer(void)
{	if(ghost==NULL && real==NULL)
		throw "double NULL ghost node";
	return ghost!=NULL ? ghost : real ;
}

