/********************************************************************************
    ElementBase.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 24 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"

// class statics
double ElementBase::gridTolerance=1.e20;

// Element globals
ElementBase **theElements;	// list of elements
int nelems=0;			// number of elements

#ifdef MPM_CODE
int ElementBase::useGimp=POINT_GIMP;
int ElementBase::gridNiNodes=4;
int ElementBase::numCPDINodes=4;
double ElementBase::rcrit=-1.;
#endif

#pragma mark ElementBase::Methods

//	Print number, ID, and all nodes to output MPM file
void ElementBase::PrintElement(ostream &os,int np)
{
    char eline[50];
	size_t esize=50;
    int i;
    int numnds=NumberNodes();
    int elemID=ElementName();
    
    // write num and elemID
#ifdef MPM_CODE
	snprintf(eline,esize,"%5d %2d  ",num,elemID);
#else
    if(np==AXI_SYM)
    {   snprintf(eline,esize,"%5d %2d %2d  %7.2lf            ",
                            num,elemID,material,GetAngleInDegrees());
    }
    else  
    {	snprintf(eline,esize,"%5d %2d %2d  %7.2lf %10.5lf ",
                            num,elemID,material,GetAngleInDegrees(),GetThickness());
    }
#endif
    os << eline;
    
    // write each node number
    for(i=0;i<numnds;i++)
    {	snprintf(eline,esize,"%5d ",nodes[i]);
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
void ElementBase::FindCentroid(Vector *center) const
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
double ElementBase::GetThickness(void) const { return (double)1.; }
void ElementBase::SetThickness(double thick) { }

// centroid (but possibly may not be in the element if it is distorted
void ElementBase::GetXYZCentroid(Vector *center)
{	center->x=(xmin+xmax)/2.;
	center->y=(ymin+ymax)/2.;
	center->z=0.;
}

// depth - 3D element return z extent, 2D element thickness in MPM not set so defaults to 1
Vector ElementBase::GetDeltaBox(void) const { return MakeVector(xmax-xmin,ymax-ymin,GetThickness()); }
double ElementBase::GetCenterX(void) const { return 0.5*(xmax+xmin); }
double ElementBase::GetCenterY(void) const { return 0.5*(ymax+ymin); }
double ElementBase::GetDeltaX(void) const { return xmax-xmin; }
double ElementBase::GetDeltaY(void) const { return ymax-ymin; }
double ElementBase::GetDeltaZ(void) const { return GetThickness(); }
bool ElementBase::IntersectsBox(Vector orig,double xlength,double ylength) const
{   if(xmax<orig.x) return false;
	if(xmin>orig.x+xlength) return false;
	if(ymax<orig.y) return false;
	if(ymin>orig.y+ylength) return false;
	return true;
}

// 2D range (3D must override)
void ElementBase::GetRange(int ax,double &amin,double &amax) const
{	if(ax==0)
	{	amin = xmin;
		amax = xmax;
	}
	else
	{	amin = ymin;
		amax = ymax;
	}
}

// number of sides in this element (override if differs)
int ElementBase::NumberSides(void) const { return 4; }

// return 0-based node number for node i where i is 1 to numnds
int ElementBase::NodeIndex(int i) { return nodes[i-1]-1; }

#pragma mark CLASS METHODS

// Return length of the shortest element side (because gridTolerance = TOLERANCE_RATIO*(shortest side)
double ElementBase::GetMinimumCellSize(void) { return gridTolerance/TOLERANCE_RATIO; }

#ifdef PREHASH_CRACKS
// Add CrackNum to crack list on this element. If no
// list yeat, create it and add this crack. If has last
// make sure CrackNum is not same as previous one.
// This method assumes cracks processed sequentially
void ElementBase::PushCrackNumOnList(int CrackNum)
{   if(SeesCrack == nullptr)
    {   SeesCrack = new vector<int>;
		SeesCrack->push_back(CrackNum);
	}
	else if(CrackNum != SeesCrack->back())
        SeesCrack->push_back(CrackNum);
}

// For next time step, delete crack list if present
void ElementBase::DeleteCrackList(void)
{   if(SeesCrack != nullptr)
    {   delete SeesCrack;
		SeesCrack = nullptr;
	}
}
#endif

#ifdef MPM_CODE
#define MAX_NODE_CONNECTIONS 10
#else
#define MAX_NODE_CONNECTIONS 41
#endif

typedef struct {
	int degree;
	int cons[MAX_NODE_CONNECTIONS];
} Connections;

// Generate mappings to resequence the nodes start at node resequence (1 based)
// On return:
//		(*outNodeMap)[j] (j=0 to nnodes-1) is original node number (0 based)
//					(can by NULL if not needed)
//		(*outRevMap)[j] is reordered node number (1 based) for original node j (1 based)
// throws std::bad_alloc, if error returns a message, not error returns NULL
char *ElementBase::ReverseMapNodes(int resequence,int **outNodeMap,int **outRevMap)
{
	int i,j,k,l,numnds;
	
	// set up data structures
	Connections *nList=new Connections[nnodes];	// Nodal connectivities
	int *theLevel=new int[nnodes];				// list of nodes in a level
	int *lastLevel=new int[nnodes];				// list of nodes in previous level
	bool *levelFlags=new bool[nnodes];			// flags to remember nodes in level
	bool *mapFlags=new bool[nnodes];			// bits to remember nodes that have already been mapped
	int *nodeMap=new int[nnodes];				// map of resequenced nodes
	
	for(i=0;i<nnodes;i++)
	{	nList[i].degree=0;			// zero the degrees
		mapFlags[i]=false;			// clear all flags
	}

	/* Invert element list to nodal connectivity list (note list has node #'s-1)
		nList[0 to nnode-1] says that node is connected to .degree nodes
												 listed in .cons[0] to .degree-1
	*/
	int node1,node2,degree;
	bool addNode;
	for(i=0;i<nelems;i++)
	{	// loop over all nodes in this element
		numnds=theElements[i]->NumberNodes();
		for(j=1;j<=numnds;j++)
		{	// get zero-based node number in the element
			node1=theElements[i]->NodeIndex(j);
			degree=nList[node1].degree;
			
			// loop over nodes in this element, skipped j
			for(k=1;k<=numnds;k++)
			{	if(k==j) continue;
				
				// zero based node number add to list if not already there
				node2=theElements[i]->NodeIndex(k);
				addNode=true;
				for(l=0;l<degree;l++)
				{	if(node2==nList[node1].cons[l])
					{	addNode=false;
						break;
					}
				}
				if(!addNode) continue;
				
				// add to list
				if(degree>=MAX_NODE_CONNECTIONS)
				{	char *msg = new char[100];
					size_t msgSize=100;
					snprintf(msg,msgSize,"Mesh too highly connected for resequencing at node %d.",node1+1);
					return msg;
				}
				nList[node1].cons[degree++]=node2;
			}
			
			// update number connected to node1
			nList[node1].degree=degree;
		}
	}
	
	// variables
	int mapped,inLastLevel,inLevel,lognb2;
	
	// Start node map using requested resequence node (-1 from 1 based to zero based here)
	nodeMap[0]=resequence-1;
	mapped=1;
	mapFlags[nodeMap[0]]=true;
	
	// Initialize level structure
	inLastLevel=1;
	lastLevel[0]=nodeMap[0];
	
	// Loop through levels until map is done
	while(true)
	{	// zero the level flags
		for(i=0;i<nnodes;i++) levelFlags[i]=false;
		
		// Calculate new level from last level
		inLevel=0;
		for(i=0;i<inLastLevel;i++)
		{	// loop over nodes connect to lode in last levenb
			node1=lastLevel[i];
			for(j=0;j<nList[node1].degree;j++)
			{	node2=nList[node1].cons[j];
				if(!mapFlags[node2] && !levelFlags[node2])
				{	theLevel[inLevel++]=node2;
					levelFlags[node2]=true;
				}
			}
		}
		
		// If inLevel is zero then all done
		if(inLevel==0) break;
		
		// Sort by degree - shell sort - Numerical Recipes in C, pg 244
		lognb2=(int)(log((double)inLevel)*1.442695022+1.0e-5);	// log base 2
		k=inLevel;
		for(l=1;l<=lognb2;l++)
		{	k>>=1;		// divide by 2
			for(j=k;j<inLevel;j++)
			{	i=j-k;
				degree=nList[theLevel[j]].degree;
				node1=theLevel[j];
				
				// Back up until find insertion point
				while(i>=0 && nList[theLevel[i]].degree>degree)
				{	theLevel[i+k]=theLevel[i];
					i-=k;
				}
				
				// Insert point
				theLevel[i+k]=node1;
			}
		}
		
		// Assign maps to nodes
		for(i=0;i<inLevel;i++)
		{	nodeMap[mapped++]=theLevel[i];
			mapFlags[theLevel[i]]=true;
		}
		
		// Copy level to last level
		inLastLevel=inLevel;
		for(i=0;i<inLevel;i++)
			lastLevel[i]=theLevel[i];
	}
	
	// delete data structures (except nodeMap)
	delete [] nList;
	delete [] theLevel;
	delete [] lastLevel;
	delete [] mapFlags;
	delete [] levelFlags;
	
	// verify found all nodes
	if(mapped!=nnodes)
	{	char *msg = new char[200];
		strcpy(msg,"Resequencing failed because of disconnected nodes. Turn resequencing off and plot the mesh to find them,\n  watching out for areas that touch but do not share nodes.");
		return msg;
	}
	
	// reverse the node map
	inLevel=(nnodes-2)/2;
	for(i=0;i<=inLevel;i++)
	{	k=nodeMap[i];
		j=nnodes-1-i;
		nodeMap[i]=nodeMap[j];
		nodeMap[j]=k;
	}
	
	// Now nodeMap[j] for j from 0 to nnodes-1 means node j (zero based) in
	// reordered nodes is nodeMap[j] in the initial listing
	
	// Calculate reverse mapping such that given node j (1 based) in the original
	// numbering, revMap[j] is node number (1 based) in the reordered list
	int *revMap=new int[nnodes+1];
	for(i=1;i<=nnodes;i++)
		revMap[nodeMap[i-1]+1]=i;
	
	// return forward and reverse mappings
	if(outNodeMap==NULL)
		delete [] nodeMap;
	else
		*outNodeMap = nodeMap;
	*outRevMap = revMap;
	
	// no errors
	return NULL;
}


