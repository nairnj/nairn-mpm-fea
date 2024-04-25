/********************************************************************************
	GridNode.cpp
	nairn-mpm-fea

	Created by John Nairn on 4/11/13.
	Copyright (c) 2013 John A. Nairn, All rights reserved.
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Patches/GhostNode.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/TransportTask.hpp"
#include "Exceptions/CommonException.hpp"

#pragma mark INITIALIZATION

// Constructors
// Create for node at zero-based (row,col) grid position
// throws std::bad_alloc
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
		ghost = NodalPoint::CreateGhostFromReal(real);
	else
		ghost = NULL;
		
	//cout << "...... node = " << realNode << "," << ghost << "," << real << "," << (real==ghost) << endl;
}

// Create for node at zero-based (row,col) grid position
// throws std::bad_alloc
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
		ghost = NodalPoint::CreateGhostFromReal(real);
	else
		ghost = NULL;
	
	//cout << "...... node = " << realNode << "," << ghost << "," << real << "," << (real==ghost) << endl;
}


#pragma mark GhostNode: Methods

// initialize ghost nodes for next time step
void GhostNode::InitializeForTimeStep(void)
{	if(ghost!=NULL)
		ghost->InitializeForTimeStep();
}

// copy allocated velocity fields from real node to ghost node
void GhostNode::InitializationReduction(void)
{	if(ghost!=NULL)
		real->CopyFieldInitialization(ghost);
}

// When mass and momentum task is done transfer ghost node values to real nodes
void GhostNode::MassAndMomentumReduction(void)
{	if(ghost!=NULL)
	{	ghost->CopyMassAndMomentum(real);
		
		// transport calculations
		TransportTask *nextTransport=transportTasks;
		while(nextTransport!=NULL)
			nextTransport=nextTransport->Task1Reduction(real,ghost);
	}
}

// Support XPIC calculations (always xpicCalculation==COPY_VSTARNEXT, and k=pass in XPIC loop)
void GhostNode::XPICSupport(int xpicCalculation,int xpicOption,NodalPoint *dummy,double timestep,int m,int k,double vsign)
{	if(ghost!=NULL)
	ghost->XPICSupport(xpicCalculation,xpicOption,real,timestep,m,k,vsign);
}

// zero gTnext on ghost nodes
void GhostNode::InitializeForXPICTransport(TransportTask *nextTransport,int xpicOption)
{	if(ghost!=NULL)
		nextTransport->InitializeForXPIC(ghost,-1.,xpicOption);
}

// XPIC reduction
void GhostNode::XPICReductionTransport(int k,TransportTask *nextTransport)
{	if(ghost!=NULL)
	{	// add from ghost to real
		TransportField *gTransGhost = nextTransport->GetTransportFieldPtr(ghost);
		TransportField *gTrans = nextTransport->GetTransportFieldPtr(real);
		gTrans->gTnext += gTransGhost->gTnext;
	
		// zero the ghost node
		gTransGhost->gTnext = 0.;
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
			nextTransport=nextTransport->ForcesReduction(real,ghost);
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
// throws CommonException()
NodalPoint *GhostNode::GetNodePointer(void)
{	if(ghost==NULL && real==NULL)
		throw CommonException("double NULL ghost node","GhostNode::GetNodePointer");
	return ghost!=NULL ? ghost : real ;
}

// return private variables as is
NodalPoint *GhostNode::GetGhostNodePointer(void) { return ghost; }
NodalPoint *GhostNode::GetRealNodePointer(void) { return real; }
