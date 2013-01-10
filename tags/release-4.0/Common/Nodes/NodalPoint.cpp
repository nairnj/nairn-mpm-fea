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

/********************************************************************************
	NodalPoint: Constructors and Destructor
********************************************************************************/

NodalPoint::NodalPoint(int nodeNum)
{
    num=nodeNum;
#ifdef MPM_CODE 
	fixedDirection=0;
	cvf=NULL;
#endif
}

// Destructor
#ifdef FEA_CODE
NodalPoint::~NodalPoint()
{
}
#endif


