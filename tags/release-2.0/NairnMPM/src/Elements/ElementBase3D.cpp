/********************************************************************************
    ElementBase3D.cpp
    NairnMPM
    
    Created by John Nairn on 7/20/06.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Elements/ElementBase3D.hpp"
#include "Nodes/NodalPoint.hpp"
#ifdef MPM_CODE
	#include "NairnMPM_Class/MeshInfo.hpp"
#endif

#pragma mark ElementBase3D::Constructors and Destructors

// main constructor - pass up the chain
ElementBase3D::ElementBase3D(long eNum,long *eNode) : ElementBase(eNum,eNode)
{
}

#pragma mark ElementBase3D::Methods

//	Find extent of this element - called once at start (and must be called)
void ElementBase3D::FindExtent(void)
{
	// find x-y extent in super class
	ElementBase::FindExtent();
	
    // find z extent of element
    int i,numnds=NumberNodes();
	double zNode;
    zmin=zmax=nd[nodes[0]]->z;
    for(i=1;i<numnds;i++)
    {   zNode=nd[nodes[i]]->z;
        if(zNode>zmax) zmax=zNode;
        if(zNode<zmin) zmin=zNode;
    }
	
    // set grid tolerance (1/10 minimum grid spacing)
    double range=TOLERANCE_RATIO*(zmax-zmin);
    if(range<gridTolerance) gridTolerance=range;
}

#pragma mark ElementBase3D::Accessors

// face nodes not meaningful
int ElementBase3D::FaceNodes(void) { return 0; }

// centroid (but possibly may not be in the element if it is distorted
void ElementBase3D::GetXYZCentroid(Vector *center)
{	center->x=(xmin+xmax)/2.;
	center->y=(ymin+ymax)/2.;
	center->z=(zmin+zmax)/2.;
}

// depth - 3D element return z extent
double ElementBase3D::GetDeltaZ(void) { return zmax-zmin; }
bool ElementBase3D::IntersectsBox(double xorig,double yorig,double xlength,double ylength,double zslice)
{	if(zslice<zmin) return false;
	if(zslice>zmax) return false;
	return ElementBase::IntersectsBox(xorig,yorig,xlength,ylength,zslice);
}

// check if this GIMP element is on the edge of the grid
// assumes a generated 3D structured grid
bool ElementBase3D::OnTheEdge(void)
{	if(!useGimp) return FALSE;
	return mpmgrid.EdgeElement3D(num);
}

// for structured grid, return 0-terminated list of neighbors
void ElementBase3D::GetListOfNeighbors(int *theList) { mpmgrid.ListOfNeighbors3D(num,theList); }


