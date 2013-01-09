/********************************************************************************
    ElementBase.cpp
    NairnMPM
    
    Created by John Nairn on Wed Jan 24 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"

// class statics
double ElementBase::gridTolerance=1.e20;

// Element globals
ElementBase **theElements;	// list of elements
int nelems=0;			// number of elements

#ifdef MPM_CODE
int ElementBase::useGimp=POINT_GIMP;
int ElementBase::analysisGimp=POINT_GIMP;
int ElementBase::numCPDINodes=4;
#endif

#pragma mark ElementBase::Methods

//	Print number, ID, and all nodes to output MPM file
void ElementBase::PrintElement(ostream &os,int np)
{
    char eline[50];
    int i;
    int numnds=NumberNodes();
    int elemID=ElementName();
    
    // write num and elemID
#ifdef MPM_CODE
	sprintf(eline,"%5d %2d  ",num,elemID);
#else
    if(np==AXI_SYM)
    {   sprintf(eline,"%5d %2d %2d  %7.2lf            ",
                            num,elemID,material,GetAngleInDegrees());
    }
    else  
    {	sprintf(eline,"%5d %2d %2d  %7.2lf %10.5lf ",
                            num,elemID,material,GetAngleInDegrees(),GetThickness());
    }
#endif
    os << eline;
    
    // write each node number
    for(i=0;i<numnds;i++)
    {	sprintf(eline,"%5d ",nodes[i]);
        os << eline;
    }
    os << endl;
}

//	Find extent of this element - called once at start (and must be called)
void ElementBase::FindExtent(void)
{	
    int i,numnds=NumberNodes();
    double xNode,range;
    
    // find extent of element
    xmin=xmax=nd[nodes[0]]->x;
    ymin=ymax=nd[nodes[0]]->y;
    for(i=1;i<numnds;i++)
	{	xNode=nd[nodes[i]]->x;
        if(xNode>xmax) xmax=xNode;
        if(xNode<xmin) xmin=xNode;
        xNode=nd[nodes[i]]->y;
        if(xNode>ymax) ymax=xNode;
        if(xNode<ymin) ymin=xNode;
    }
	
    // set grid tolerance (1/10 minimum grid spacing)
    range=TOLERANCE_RATIO*(xmax-xmin);
    if(range<gridTolerance) gridTolerance=range;
    range=TOLERANCE_RATIO*(ymax-ymin);
    if(range<gridTolerance) gridTolerance=range;
}

//	Find center of mass of element (2D), and needed before extent is known
void ElementBase::FindCentroid(Vector *center)
{
    int i,numnds=NumberNodes();
	double xtot=nd[nodes[0]]->x;
	double ytot=nd[nodes[0]]->y;
    for(i=1;i<numnds;i++)
	{	xtot+=nd[nodes[i]]->x;
        ytot+=nd[nodes[i]]->y;
    }
	center->x=xtot/(double)numnds;
	center->y=ytot/(double)numnds;
	center->z=0.;
}

#pragma mark ElementBase::Accessors

// thickness which may be in a subclass
double ElementBase::GetThickness(void) { return (double)1.; }
void ElementBase::SetThickness(double thick) { }

// centroid (but possibly may not be in the element if it is distorted
void ElementBase::GetXYZCentroid(Vector *center)
{	center->x=(xmin+xmax)/2.;
	center->y=(ymin+ymax)/2.;
	center->z=0.;
}

// depth - 3D element return z extent
double ElementBase::GetCenterX(void) { return 0.5*(xmax+xmin); }
double ElementBase::GetDeltaX(void) { return xmax-xmin; }
double ElementBase::GetDeltaY(void) { return ymax-ymin; }
double ElementBase::GetDeltaZ(void) { return GetThickness(); }
bool ElementBase::IntersectsBox(double xorig,double yorig,double xlength,double ylength,double zslice)
{	if(xmax<xorig) return false;
	if(xmin>xorig+xlength) return false;
	if(ymax<yorig) return false;
	if(ymin>yorig+ylength) return false;
	return true;
}

// number of sides in this element (override if differs)
int ElementBase::NumberSides(void) { return 4; }

// return 0-based node number for node i where i is 1 to numnds
int ElementBase::NodeIndex(int i) { return nodes[i-1]-1; }

#pragma mark CLASS METHODS

// Return length of the sortest element side (because gridTolerance = TOLERANCE_RATIO*(shortest side)
double ElementBase::GetMinimumCellSize(void) { return gridTolerance/TOLERANCE_RATIO; }



