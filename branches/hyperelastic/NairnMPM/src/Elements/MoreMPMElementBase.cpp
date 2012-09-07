/********************************************************************************
    MoreMPMElementBase.cpp
    NairnMPM
    
    Created by John Nairn on Fri Oct 22 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"

#define MAXITER 100


#pragma mark ElementBase: Constructors and Destructor MPM Only

/* Main MPM contructor when creating elements:
	input is node numbers (0-based array) but values
            are 1-based (length always MaxElNd, even if unused)
        subclass should set nodes[NumberNodes]=nodes[0]
*/
ElementBase::ElementBase(int eNum,int *eNode)
{	int i;
    num=eNum;
    for(i=0;i<MaxElNd;i++) nodes[i]=eNode[i];
	neighbors=NULL;
    filled=0;
}

ElementBase::~ElementBase()
{	if(neighbors!=NULL) delete [] neighbors;
}

#pragma mark ElementBase: ShapeFunction calls for MPM

/* Find dimensionless position first, then get shape functions
	Load number of nodes into numnds
	Load dimensionless position into xipos vector
	Load node numbers into nds[1]...
	Load shape functions into fn[1]...
	Input: pointer to material point position and dimensionless position
	See other GetShapeFunctions() if need to change
*/
void ElementBase::GetShapeFunctions(int *numnds,double *fn,int *nds,Vector *pos,Vector *xipos)
{
	if(useGimp)
	{	// GIMP analysis
		int ndIDs[MaxShapeNds];
		
		GetXiPos(pos,xipos);
		GetGimpNodes(numnds,nds,ndIDs,xipos);
		GimpShapeFunction(xipos,*numnds,ndIDs,FALSE,&fn[1],NULL,NULL,NULL);
		GimpCompact(numnds,nds,fn,NULL,NULL,NULL);
	}
	else
	{	// Load element noodes, dimensionless position, and shape functinos
		GetNodes(numnds,nds);
		GetXiPos(pos,xipos);
		ShapeFunction(xipos,FALSE,&fn[1],NULL,NULL,NULL);
	}
}

/* Just get nodes and shape functions (when know dimensionless position)
	Load number of nodes into numnds
	Load node numbers into nds[1]...
	Load shape functions into fn[1]...
	Input: pointer to material point dimensionless position
	See other GetShapeFunctions() if need to change
*/
void ElementBase::GetShapeFunctions(int *numnds,double *fn,int *nds,Vector *xipos)
{
	if(useGimp)
	{	// GIMP analysis
		int ndIDs[MaxShapeNds];
		
		GetGimpNodes(numnds,nds,ndIDs,xipos);
		GimpShapeFunction(xipos,*numnds,ndIDs,FALSE,&fn[1],NULL,NULL,NULL);
		GimpCompact(numnds,nds,fn,NULL,NULL,NULL);
	}
	else
	{	// load coordinates if not already done
		GetNodes(numnds,nds);
		ShapeFunction(xipos,FALSE,&fn[1],NULL,NULL,NULL);
	}
}

/* do several element things at once
	Load number of nodes into numnds
	Load node numbers into nds[1]...
	Load shape functions into fn[1]...
	Load shape function derviatives into xDeriv[1]..., yDeriv[1]..., zDeriv[1]...
		zDeriv may be NULL for 2D
	Input: pointer to material point dimensionless position
	See all GetShapeFunctionsAndGradients() if need to change
*/
void ElementBase::GetShapeGradients(int *numnds,double *fn,int *nds,Vector *xipos,
										double *xDeriv,double *yDeriv,double *zDeriv)
{
	if(useGimp)
	{	// GIMP analysis
		int ndIDs[MaxShapeNds];
		
		GetGimpNodes(numnds,nds,ndIDs,xipos);
		GimpShapeFunction(xipos,*numnds,ndIDs,TRUE,&fn[1],&xDeriv[1],&yDeriv[1],&zDeriv[1]);
		GimpCompact(numnds,nds,fn,xDeriv,yDeriv,zDeriv);
	}
	else
	{	// load nodal numbers
		GetNodes(numnds,nds);
	
		// special case for regular mesh
		if(mpmgrid.GetCartesian()>0)
		{	double *zptr = zDeriv==NULL ? NULL : &zDeriv[1];
			ShapeFunction(xipos,TRUE,&fn[1],&xDeriv[1],&yDeriv[1],zptr);
		}
		else
		{	// Load element coordinates
			Vector ce[MaxElNd];
			double fnh[MaxElNd];
			GetCoordinates(ce,*numnds,nds);

			// find shape functions and derviatives
			ShapeFunction(xipos,BMATRIX,&fn[1],&xDeriv[1],&yDeriv[1],&ce[1],NULL,NULL,&fnh[1]);
		}
	}
}

/* Do several element things at once including finding dimensionless position of current particle
	only called in MPM and multimaterial mode at start of time step
	Load number of nodes into numnds
	Load dimensionless position of particle into xipos on the particle
	Load node numbers into nds[1]...
	Load shape functions into fn[1]...
	Load shape function derviatives into xDeriv[1]..., yDeriv[1]..., zDeriv[1]...
	zDeriv may be NULL for 2D
	Input: pointer to material point dimensionless position, which is changed here
	See also GetShapeGradients() if need to change
*/
void ElementBase::GetShapeFunctionsAndGradients(int *numnds,double *fn,int *nds,Vector *pos,Vector *xipos,
									double *xDeriv,double *yDeriv,double *zDeriv)
{
	if(useGimp)
	{	// GIMP analysis
		int ndIDs[MaxShapeNds];
		
		GetXiPos(pos,xipos);
		GetGimpNodes(numnds,nds,ndIDs,xipos);
		GimpShapeFunction(xipos,*numnds,ndIDs,TRUE,&fn[1],&xDeriv[1],&yDeriv[1],&zDeriv[1]);
		GimpCompact(numnds,nds,fn,xDeriv,yDeriv,zDeriv);
	}
	else
	{	// load nodal numbers
		GetNodes(numnds,nds);
		GetXiPos(pos,xipos);
		
		// special case for regular mesh
		if(mpmgrid.GetCartesian()>0)
		{	double *zptr = zDeriv==NULL ? NULL : &zDeriv[1];
			ShapeFunction(xipos,TRUE,&fn[1],&xDeriv[1],&yDeriv[1],zptr);
		}
		else
		{	// Load element coordinates
			Vector ce[MaxElNd];
			double fnh[MaxElNd];
			GetCoordinates(ce,*numnds,nds);
			
			// find shape functions and derviatives
			ShapeFunction(xipos,BMATRIX,&fn[1],&xDeriv[1],&yDeriv[1],&ce[1],NULL,NULL,&fnh[1]);
		}
	}
}

/* Find dimensionless coordinates by numerical methods
	input: pos is position in the element
	output: xipos is dimensionless position
	only used in MPM and only here if non-rectangular elements
*/
void ElementBase::GetXiPos(Vector *pos,Vector *xipos)
{
    double xt,yt,dxxi,dxeta,dyxi,dyeta;
    double deter,dxi,deta,dist;
    double gfn[MaxElNd],gdfnxi[MaxElNd],gdfnet[MaxElNd];
    int numnds=NumberNodes(),i,j;
	
	// initial guess
	GetCentroid(xipos);
	
	// nodal coordinates
	Vector eNode[MaxElNd];
	for(i=0;i<numnds;i++)
	{	eNode[i].x=nd[nodes[i]]->x;
		eNode[i].y=nd[nodes[i]]->y;
	}
    
    /* solve for xipos using Newton-Rapheson (see FEA Notes)
            using shape functions and their derivatives */
    for(j=1;j<=MAXITER;j++)
    {   ShapeFunction(xipos,TRUE,gfn,gdfnxi,gdfnet,NULL,NULL,NULL,NULL);
        xt=-pos->x;
        yt=-pos->y;
        dxxi=0.;
        dxeta=0.;
        dyxi=0.;
        dyeta=0.;
        for(i=0;i<numnds;i++)
        {   xt+=eNode[i].x*gfn[i];
            yt+=eNode[i].y*gfn[i];
            dxxi+=eNode[i].x*gdfnxi[i];
            dxeta+=eNode[i].x*gdfnet[i];
            dyxi+=eNode[i].y*gdfnxi[i];
            dyeta+=eNode[i].y*gdfnet[i];
        }
        deter=dxxi*dyeta-dxeta*dyxi;
        dxi=(-xt*dyeta+yt*dxeta)/deter;
        deta=(xt*dyxi-yt*dxxi)/deter;
        xipos->x+=dxi;
        xipos->y+=deta;
        dist=sqrt(dxi*dxi+deta*deta);
        if(dist<.001) break;
    }
}

// return dimensionless location for material points
void ElementBase::MPMPoints(short numPerElement,Vector *mpos)
{
    int j,k;
    double fxn[MaxElNd];
    
    if(NumberSides()==4)
    {	switch(numPerElement)
        {   case 4:
                // ENI or FNI - 2D only
                mpos[0].x=-.5;
                mpos[0].y=-.5;
                mpos[0].z=0.;
                mpos[1].x=.5;
                mpos[1].y=-.5;
                mpos[1].z=0.;
                mpos[2].x=-.5;
                mpos[2].y=.5;
                mpos[2].z=0.;
                mpos[3].x=.5;
                mpos[3].y=.5;
                mpos[3].z=0.;
                break;
            case 1:
                // CM of square or brick
                mpos[0].x=0.;
                mpos[0].y=0.;
				mpos[0].z=0.;
				break;
			case 8:
				// 3D box
                mpos[0].x=-.5;
                mpos[0].y=-.5;
				mpos[0].z=-.5;
                mpos[1].x=.5;
                mpos[1].y=-.5;
				mpos[1].z=-.5;
                mpos[2].x=-.5;
                mpos[2].y=.5;
				mpos[2].z=-.5;
                mpos[3].x=.5;
                mpos[3].y=.5;
				mpos[3].z=-.5;
                mpos[4].x=-.5;
                mpos[4].y=-.5;
				mpos[4].z=.5;
                mpos[5].x=.5;
                mpos[5].y=-.5;
				mpos[5].z=.5;
                mpos[6].x=-.5;
                mpos[6].y=.5;
				mpos[6].z=.5;
                mpos[7].x=.5;
                mpos[7].y=.5;
				mpos[7].z=.5;
				break;
            default:
                break;
        }
    }
    
    // covert to x-y-z locations
    for(k=0;k<numPerElement;k++)
    {	ShapeFunction(&mpos[k],FALSE,fxn,NULL,NULL,NULL);
		ZeroVector(&mpos[k]);
        for(j=0;j<NumberNodes();j++)
        {   mpos[k].x+=nd[nodes[j]]->x*fxn[j];
            mpos[k].y+=nd[nodes[j]]->y*fxn[j];
            mpos[k].z+=nd[nodes[j]]->z*fxn[j];
        }
    }
}

// Get GIMP nodes around an element #num, but only where shape functions is non zero
// assumed to be properly numbered regular array
// load nodes into nds[1]... and node IDs (0-15) into ndIDs[0]...
// Elements that support GIMP must override
void ElementBase::GetGimpNodes(int *numnds,int *nds,int *ndIDs,Vector *xipos)
{
}

// get GIMP shape functions and optionally derivatives wrt x and y
// assumed to be properly numbered regular array
// input *xi position in element coordinate and ndIDs[0]... is which nodes (0-15)
// Elements that support GIMP must override
void ElementBase::GimpShapeFunction(Vector *xi,int numnds,int *ndIDs,int getDeriv,double *sfxn,
										 double *xDeriv,double *yDeriv,double *zDeriv)
{
}

// Ignore zero shape functions in GIMP, but only need to check remote nodes
void ElementBase::GimpCompact(int *numnds,int *nds,double *sfxn,double *xDeriv,double *yDeriv,double *zDeriv)
{
	int local=NumberNodes();
	int i,found=local;
	for(i=local+1;i<=*numnds;i++)
	{	if(sfxn[i]>0.0)
		{	found++;
			nds[found]=nds[i];
			sfxn[found]=sfxn[i];
			if(xDeriv!=NULL)
			{	xDeriv[found]=xDeriv[i];
				yDeriv[found]=yDeriv[i];
				if(zDeriv!=NULL) zDeriv[found]=zDeriv[i];
			}
		}
	}
	*numnds=found;
}

// check if this GIMP element is on the edge of the grid
// assumes a generated 2D structured grid
bool ElementBase::OnTheEdge(void)
{	if(!useGimp) return FALSE;
	return mpmgrid.EdgeElement2D(num);
}

#pragma mark ElementBase: MPM Only Methods

/* Get list of nodes for this element and the numbner of nodes
	nodes will be in nds[1]...
*/
void ElementBase::GetNodes(int *numnds,int *nds)
{
	int i;
	*numnds=NumberNodes();
	for(i=1;i<=*numnds;i++)
	   nds[i]=nodes[i-1];
}

/* Load nodal coordinates into eNodes[1]...
	input numnds and nds from previous call to GetNodes()
*/
void ElementBase::GetCoordinates(Vector *eNodes,int numnds,int *nds)
{
	int i;
	for(i=1;i<=numnds;i++)
	{	eNodes[i].x=nd[nds[i]]->x;
		eNodes[i].y=nd[nds[i]]->y;
		eNodes[i].z=nd[nds[i]]->z;
	}
}

/*  Give dimensionless coordinates for central point
    of element to provide initial guess for
    numerical solution to find natural coordinates.
	Only needed for non-rectangular elements
*/
void ElementBase::GetCentroid(Vector *xipos)
{
    xipos->x=0.;
    xipos->y=0.;
    xipos->z=0.;
}

/*  Find node of this element closest to specified point
    Only checks nodes of this element. In distorted mesh,
    could be closer node in neighboring element, even if
    this point is in this element
*/
int ElementBase::NearestNode(double x,double y,int *secondChoice)
{	int ind,nearest,nextNearest;
    double dist,distMin,distNextMin,dltx,dlty;
    int k,numnds=NumberNodes();
    
    // distance to first node
	nearest=nodes[0];
    dltx=nd[nearest]->x-x;
    dlty=nd[nearest]->y-y;
    distMin=dltx*dltx+dlty*dlty;
    
    // distance to second node
	nextNearest=nodes[1];
	dltx=nd[nextNearest]->x-x;
	dlty=nd[nextNearest]->y-y;
	distNextMin=dltx*dltx+dlty*dlty;
	if(distNextMin<distMin)
	{	dist=distMin;
		distMin=distNextMin;
		distNextMin=distMin;
		nearest=nodes[1];
		nextNearest=nodes[0];
	}
	
    // check the rest
    for(k=2;k<numnds;k++)
    {	ind=nodes[k];
    	dltx=nd[ind]->x-x;
        dlty=nd[ind]->y-y;
        dist=dltx*dltx+dlty*dlty;
        if(dist<distMin)
		{	distNextMin=distMin;
			nextNearest=nearest;
			distMin=dist;
            nearest=ind;
        }
		else if(dist<distNextMin)
		{	distNextMin=dist;
			nextNearest=ind;
		}
    }
    
    // return the node
	*secondChoice=nextNearest;
    return nearest;
}

/* find next node (in ccw direction) from specified node
    Assumes nodes[NumberNodes()]=nodes[0]
    return 0 if gridNode not in this element
*/
int ElementBase::NextNode(int gridNode)
{	int i;
    for(i=0;i<NumberNodes();i++)
    {	if(nodes[i]==gridNode)
            return nodes[i+1];
    }
    return 0;
}

// for structured grid, return 0-terminated list of neighbors
void ElementBase::GetListOfNeighbors(int *theList) { mpmgrid.ListOfNeighbors2D(num,theList); }

// If needed, create the neighbors array (only down for 2D with cracks)
void ElementBase::AllocateNeighborsArray(void)
{	neighbors = new int[NumberSides()];
	int i;
	for(i=0;i<NumberSides();i++)
		neighbors[i]=UNKNOWN_NEIGHBOR;
}

/* Find element (0 based) that touches the element face stating in gridNode
    First time called, search other elements for the edge, but then store result
    Assumes nodes[NumberNodes()]=nodes[0]
    return NO_NEIGHBOR (-1) if no neighbor (i.e., on edge of grid)
    return UNKNOWN_NEIGHBOR (-2) if gridNode not in this element
*/
int ElementBase::Neighbor(int gridNode)
{
    int i,edge=-1;
    
    // which edge?
    for(i=0;i<NumberSides();i++)
    {	if(nodes[i]==gridNode)
        {   edge=i;
            break;
        }
    }
    if(edge<0) return UNKNOWN_NEIGHBOR;
    
    // search for neighbor if needed
    if(neighbors[edge]==UNKNOWN_NEIGHBOR)
    {	int nextNode=nodes[edge+1];
        neighbors[edge]=NO_NEIGHBOR;
        
        // search for element with nextNode,gridNode in ccw direction
        for(i=0;i<nelems;i++)
        {   if(theElements[i]->FindEdge(nextNode,gridNode)>=0)
            {	neighbors[edge]=i;
                break;
            }
        }
    }
    
    // return result
    return neighbors[edge];
}

/* find an edge (in ccw direction) with a pair of nodes
    Assumes nodes[NumberNodes()]=nodes[0]
    return edge number (0 based) or -1 if no such edge
*/
int ElementBase::FindEdge(int beginNode,int endNode)
{
    int i;
    for(i=0;i<NumberNodes();i++)
    {	if(nodes[i]==beginNode && nodes[i+1]==endNode)
            return i;
    }
    return -1;
}

/* return TRUE or false if Orthogonal in x-y-z planes and if orthogonal then find width, height, and depth
	of the element. If 2D element, return depth=0
	*/
int ElementBase::Orthogonal(double *dx,double *dy,double *dz) { return FALSE; }

#pragma mark CLASS METHODS

// zero all velocity fields at start of time step
void ElementBase::AllocateNeighbors(void)
{	int i;
    for(i=0;i<nelems;i++)
        theElements[i]->AllocateNeighborsArray();
}

// turn GIMP off and then back on to get normal shape functions
// can call with out GIMP and will be no changes
void ElementBase::DisableGimp(void)
{	previousGimp=useGimp;
	useGimp=FALSE;
}
void ElementBase::EnableGimp(void) { useGimp=previousGimp; }




