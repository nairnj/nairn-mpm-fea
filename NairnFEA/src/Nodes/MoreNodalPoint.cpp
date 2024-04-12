/********************************************************************************
    MoreNodalPoint.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Mar 15 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Nodes/NodalPoint.hpp"

// allocate and initialize force field
// throws std::bad_alloc
void NodalPoint::InitForceField(void)
{
    fs=new ForceField;
    
    fs->stress.xx=0.;
    fs->stress.yy=0.;
    fs->stress.zz=0.;
    fs->stress.xy=0.;
    fs->force.x=0.;
    fs->force.y=0.;
    fs->numElems=0;
}

// allocate and initialize force field
void NodalPoint::PrintAvgStress(void)
{
    char fline[200];
    size_t fsize=200;
    double factor=1./(double)fs->numElems;
    
    snprintf(fline,fsize,"%5d  %15.7e  %15.7e  %15.7e  %15.7e",num,
            fs->stress.xx*factor,fs->stress.yy*factor,
            fs->stress.zz*factor,fs->stress.xy*factor);
    cout << fline << endl;
}

#pragma mark CLASS METHODS

// create 2D node with or without cracks
NodalPoint *NodalPoint::Create2DNode(int nodeNum,double xPt,double yPt)
{
	NodalPoint *newNode = new NodalPoint(nodeNum,xPt,yPt);
	return newNode;
}

// create 2D node with or without cracks
NodalPoint *NodalPoint::Create3DNode(int nodeNum,double xPt,double yPt,double zPt)
{
	NodalPoint *newNode = new NodalPoint(nodeNum,xPt,yPt,zPt);
	return newNode;
}
