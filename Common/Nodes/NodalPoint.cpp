/********************************************************************************
    NodalPoint.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 24 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Nodes/NodalPoint.hpp"
#ifdef MPM_CODE
	#include "NairnMPM_Class/NairnMPM.hpp"
#else
	#include "NairnFEA_Class/NairnFEA.hpp"
#endif

// Nodal Point globals
int nnodes=0;			// number of nodes
NodalPoint **nd;		// list of nodes
int *nda;				// list of active nodes (number in nda[0] = *nda

#pragma mark NodalPoint: Constructors and Destructor

// 2D nodal point
NodalPoint::NodalPoint(int nodeNum,double xPt,double yPt) : LinkedObject()
{
	CreateNodalPoint(nodeNum,xPt,yPt,0.);
}

// 3D nodal point
NodalPoint::NodalPoint(int nodeNum,double xPt,double yPt,double zPt) : LinkedObject()
{
	CreateNodalPoint(nodeNum,xPt,yPt,zPt);
}

void NodalPoint::CreateNodalPoint(int nodeNum,double xPt,double yPt,double zPt)
{
	x=xPt;
	y=yPt;
	z=zPt;
	
    num=nodeNum;
#ifdef MPM_CODE 
	fixedDirection=0;
	cvf=NULL;
	contactData=NULL;
	gDiff=NULL;
#endif
}

#ifdef MPM_CODE

// create ghost node that is a copy of source node
// throws std::bad_alloc
NodalPoint::NodalPoint(NodalPoint *real)
{
    x = real->x;
    y = real->y;
    z = real->z;
    
	num = real->num;
	fixedDirection=0;
	cvf=NULL;
    contactData=NULL;
    gDiff=NULL;

	// this called later for real nodes
	PrepareForFields();
	
	// needed for diffusion tasks
	CreateDiffusionVariables();

}

#else

// FEA Destructor
NodalPoint::~NodalPoint()
{
}

#endif

// write node to output file
void NodalPoint::PrintNodalPoint(ostream &os)
{
	char nline[200];
    size_t nlsize=200;
	if(fmobj->IsThreeD())
		snprintf(nline,nlsize,"%5d %15.7e %15.7e %15.7e",num,x,y,z);
	else
        snprintf(nline,nlsize,"%5d %15.7e %15.7e",num,x,y);
	os << nline << endl;
}

