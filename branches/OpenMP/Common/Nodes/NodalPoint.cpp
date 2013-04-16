/********************************************************************************
    NodalPoint.cpp
    NairnMPM
    
    Created by John Nairn on Wed Jan 24 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Nodes/NodalPoint.hpp"

// Nodal Point global
int nnodes=0;			// number of nodes
NodalPoint **nd;		// list of nodes

#pragma mark NodalPoint: Constructors and Destructor

NodalPoint::NodalPoint(int nodeNum)
{
    num=nodeNum;
#ifdef MPM_CODE 
	fixedDirection=0;
	cvf=NULL;
#endif
}

#ifdef MPM_CODE

// create ghost node that is a copy of source node
NodalPoint::NodalPoint(NodalPoint *real)
{
	num = real->num;
	x = real->x;
	y = real->y;
	z = real->z;
	fixedDirection=0;
	cvf=NULL;
	
	// this called later for real nodes
	PrepareForFields();
}

#else

// FEA Destructor
NodalPoint::~NodalPoint()
{
}

#endif


