/*********************************************************************
    MeshInfo.cpp
    Nairn Research Group MPM Code
    
    Created by John Nairn on 5/16/06.
    Copyright (c) 2006, All rights reserved.
********************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "NairnMPM_Class/MeshInfo.hpp"
#include "NairnMPM_Class/Reservoir.hpp"
#include "Patches/GridPatch.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Elements/ElementBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Materials/MaterialBase.hpp"
#include "Materials/ContactLaw.hpp"
#include "Exceptions/MPMWarnings.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include <algorithm>
#include <iostream>
using namespace std;
#include "System/MPMPrefix.hpp"

// global class for grid information
MeshInfo mpmgrid;
int MeshInfo::warnLRConvergence = -1;

#pragma mark MeshInfo:Constructors and Destructors

MeshInfo::MeshInfo(void)
{
	cartesian=UNKNOWN_GRID;
	horiz=0;                        // also flag that used <Grid> command (i.e. structured grid)
	cellMinSize=-1.;				// calculated when needed
	xpnum=-1;						// default patching
	tartanPatches=false;
	minParticleSize = MakeVector(1.e50,1.e50,1.e50);
	
	// for contact
	materialNormalMethod=LOGISTIC_REGRESSION;		// method to find normals in multimaterial contact
	hasImperfectInterface = false;							// flag for any imperfect interfaces
	materialContactLawID = -1;
	rigidGradientBias=1.;						// Use rigid gradient unless material volume gradient is this much higher (only normal method 2)
	lumpingMethod = LUMP_OTHER_MATERIALS;
	volumeGradientIndex = -1;		// turned on in multimaterial mode only
	positionIndex = -1;				// index to position extrapolation for contact
	displacementIndex = -1;			// index to displacement extrapolation for contact
	contactByDisplacements=true;				// contact by displacements for materials
}

#pragma mark MeshInfo:Methods

// output grid info
void MeshInfo::Output(bool isAxisym)
{
    char fline[200];
	size_t fsize=200;

	if(cartesian>0)
	{	if(horiz>0)
			cout << "Structured";
		else
			cout << "Unstructured";
		
        if(cartesian == VARIABLE_ORTHOGONAL_GRID || cartesian == VARIABLE_RECTANGULAR_GRID)
        {   cout << " orthogonal grid with variable element sizes" << endl;
        }
        
        else
        {   if(isAxisym)
				snprintf(fline,fsize," orthogonal grid with dR: %g dZ: %g",grid.x,grid.y);
            else
				snprintf(fline,fsize," orthogonal grid with dx: %g dy: %g",grid.x,grid.y);
            cout << fline;
            if(fmobj->IsThreeD())
            {	snprintf(fline,fsize," dz: %g",grid.z);
                cout << fline;
            }
            cout << endl;
            
            if(isAxisym)
				snprintf(fline,fsize,"Origin at R: %g Z: %g",xmin,ymin);
            else
				snprintf(fline,fsize,"Origin at x: %g y: %g",xmin,ymin);
            cout << fline;
            if(fmobj->IsThreeD())
            {	snprintf(fline,fsize," z: %g",zmin);
                cout << fline;
            }
            cout << endl;
        }
		
		if(!fmobj->IsThreeD() && !isAxisym)
		{	snprintf(fline,fsize,"Thickness: %g",zmin);
			cout << fline << endl;
		}

#ifdef _OPENMP
		if(isAxisym)
			snprintf(fline,fsize,"Patch Grid R: %d Z: %d",xpnum,ypnum);
		else
			snprintf(fline,fsize,"Patch Grid x: %d y: %d",xpnum,ypnum);
		cout << fline;
		if(fmobj->IsThreeD())
		{	snprintf(fline,fsize," z: %d",zpnum);
			cout << fline;
		}
		if(tartanPatches)
		{	cout << " with tartan spacings:" << endl;
			if(isAxisym)
				cout << "  R:";
			else
				cout << "  X:";
			for(int i=0;i<xsizes.size();i++)
				cout << " " << xsizes[i];
			cout << endl;
			if(isAxisym)
				cout << "  Z:";
			else
				cout << "  Y:";
			for(int i=0;i<ysizes.size();i++)
				cout << " " << ysizes[i];
			cout << endl;
			if(fmobj->IsThreeD())
			{	cout << "  Z:";
				for(int i=0;i<zsizes.size();i++)
					cout << " " << zsizes[i];
				cout << endl;
			}
		}
		else
			cout << endl;
#endif
	}
	else
		cout << "Non-orthogonal grid" << endl;
}

// output contact method by displacements or position with a cutoff
void MeshInfo::OutputContactByDisplacements(bool regressionMethod,bool byDisplacements,double cutoff)
{

	// special case for material contact when using a regression method
	if(regressionMethod)
	{	cout << "   (normal cod found in regression methods)" << endl;
		contactByDisplacements = true;		// no needed for normal, set to true for interface displacements
		return;
	}
	
	// other material contact and crack contact methods
    if(byDisplacements)
		cout << "   (normal cod from displacements)" << endl;
	else if(cutoff>0.)
	{	cout << "   (normal cod from position with contact when separated less than " << cutoff
                << " of cell)" << endl;
	}
	else
	{	cout << "   (normal cod from position with contact by separation using a " << -cutoff
						<< " power-law correction)" << endl;
	}
	
}

// check if element is on edge of 2D structured mesh - only needed for GIMP calculations
// Feature that calls this method must require the problem to have a structured <Grid>
bool MeshInfo::EdgeElement2D(int num)
{	
	// check first and last row
	if(num<=horiz || num>totalElems-horiz) return TRUE;
	
	// first or last column (remainder 0 or 1 after modulo division)
	int col=num % horiz;
	if(col<=1) return TRUE;
	
	return FALSE;
}

// check if element is on edge of 3D structured mesh - only need for GIMP calculations
// Feature that calls this method must require the problem to have a structured <Grid>
bool MeshInfo::EdgeElement3D(int num)
{
	// check x-y face
	if(num<=horiz*vert || num>totalElems-horiz*vert) return TRUE;
	
	// check y-z face (remainder 0 or 1 after modulo division)
	int yz_face=num % horiz;
	if(yz_face<=1) return TRUE;

	// check x-z face (1-horiz on front face and 0 and > h*v-h up on back face)
	int xz_face=num % (horiz*vert);
	if(xz_face<=horiz || xz_face>horiz*(vert-1)) return TRUE;
	
	return FALSE;
}

// check if node is on beyond end of row, column or rank (1-based), dir is 'x','y', 'z'
// Feature that calls this method must require the problem to have a structured <Grid>
bool MeshInfo::EdgeNode(int num,char dir)
{
	switch(dir)
	{	case 'x':
			return num >= horiz+1;
		case 'y':
			return num >= vert+1;
		case 'z':
			return num >= depth+1;
		default:
			break;
	}
	return TRUE;
}

// If have a structured grid, get the eight 2D neighbor elements. The 1-based
// element numbers are returned and the list is terminated by 0
// unless on the edge, it will return 8 (3 below, left, right, 3 above)
// neighbor needs to be size [9]
void MeshInfo::ListOfNeighbors2D(int num,int *neighbor)
{
	if(horiz<=0)
	{	neighbor[0]=0;
		return;
	}
	
	// find the column (here 0 is last column and 1 is first column)
	// element num = (row-1)*Nhoriz + col;		// 1 based
	int i=0;
	int col=num % horiz;
	int below=num-horiz;			// number in row below
	int above=num+horiz;			// number in row above
	
	if(col==1)
	{	// the element in is the first column having up to 5 neighbors
		if(num>horiz)
		{	// elements in row below
			neighbor[i++]=below;
			neighbor[i++]=below+1;
		}
		// elements in same row
		neighbor[i++]=num+1;
		if(num<=totalElems-horiz)
		{	// elements in row above
			neighbor[i++]=above+1;
			neighbor[i++]=above;
		}
	}
		
	else if(col==0)
	{	// the element in is the last column having up to 5 neighbors
		if(num>horiz)
		{	// elements in row below
			neighbor[i++]=below-1;
			neighbor[i++]=below;
		}
		// elements in same row
		neighbor[i++]=num-1;
		if(num<=totalElems-horiz)
		{	// elements in row above
			neighbor[i++]=above-1;
			neighbor[i++]=above;
		}
	}
	
	else
	{	// the element is not on the edge
		if(num>horiz)
		{	// elements in row below
			neighbor[i++]=below-1;
			neighbor[i++]=below;
			neighbor[i++]=below+1;
		}
		// elements in same row
		neighbor[i++]=num-1;
		neighbor[i++]=num+1;
		if(num<=totalElems-horiz)
		{	// elements in row above
			neighbor[i++]=above-1;
			neighbor[i++]=above;
			neighbor[i++]=above+1;
		}
	}
	
	// mark the end
	neighbor[i]=0;
}

// If have a structured, orthogonal grid, get elements -1 in all direction in ebelow
// and +1 in all direction in eabove (used by Tartan grid shape functions)
// Returns 0-based numbers to point to element list
// Warning: does not check and below>=0 and above<umber of elements
void MeshInfo::GetCorners2D(int num,int &ebelow,int &eabove)
{
    // element num = (row-1)*Nhoriz + col;        // 1 based
    ebelow = num-2-horiz;
    eabove = num+horiz;
}

// If have a structured grid, get the 26 (3D) neighbor elements. The 1-based
// element numbers are returned and the list is terminated by 0
// neighbor needs to be size [27]
void MeshInfo::ListOfNeighbors3D(int num,int *neighbor)
{	// exit if not structure grid
	if(horiz<=0)
	{	neighbor[0]=0;
		return;
	}
	
	// col is constant x, row is constant y, and slice is constant z
	// element # = (slice-1)*Nhoriz*Nvert + (row-1)*Nhoriz + col;  // if all 1 based
	int i=0;
	int perSlice = horiz*vert;			// number in each slice
	int snum = (num-1) % perSlice;      // index in element's slice (0 to horiz*vert-1)
	int col = snum % horiz;				// x column for this element (0 to horiz-1)
	
	// Do each slice
	int j;
    // 1-based element numbers below, at, and above element num
	for(j=num-perSlice;j<=num+perSlice;j+=perSlice)
	{	// element numbers in neighboring rows in this slice
		if(j<1 || j>totalElems) continue;
		int below = j-horiz;			// element # in row below or above
		int above = j+horiz;
		
		if(col==0)
		{	// the element is in the first column having up to 5 neighbors
			if(snum>=horiz)
			{	// elements in row below
				neighbor[i++]=below;
				neighbor[i++]=below+1;
			}
			// elements in same row
			if(j!=num) neighbor[i++]=j;
			neighbor[i++]=j+1;
			if(snum<perSlice-horiz)
			{	// elements in row above
				neighbor[i++]=above+1;
				neighbor[i++]=above;
			}
		}
		
		else if(col==horiz-1)
		{	// the element in is the last column having up to 5 neighbors
			if(snum>=horiz)
			{	// elements in row below
				neighbor[i++]=below-1;
				neighbor[i++]=below;
			}
			// elements in same row
			neighbor[i++]=j-1;
			if(j!=num) neighbor[i++]=j;
			if(snum<perSlice-horiz)
			{	// elements in row above
				neighbor[i++]=above-1;
				neighbor[i++]=above;
			}
		}
		
		else
		{	// the element is not on the edge
			if(snum>=horiz)
			{	// elements in row below
				neighbor[i++]=below-1;
				neighbor[i++]=below;
				neighbor[i++]=below+1;
			}
			// elements in same row
			neighbor[i++]=j-1;
			if(j!=num) neighbor[i++]=j;
			neighbor[i++]=j+1;
			if(snum<perSlice-horiz)
			{	// elements in row above
				neighbor[i++]=above-1;
				neighbor[i++]=above;
				neighbor[i++]=above+1;
			}
		}
	}
	
	// mark the end
	neighbor[i]=0;
}

// If have a structured, orthogonal grid, get elements -1 in all direction in ebelow
// and +1 in all direction in eabove (used by Tartan grid shape functions)
// Returns 0-based numbers to point to element list
// Warning: does not check and below>=0 and above<umber of elements
void MeshInfo::GetCorners3D(int num,int &ebelow,int &eabove)
{
    // element num = (slice-1)*Nhoriz*Nvert + (row-1)*Nhoriz + col;  // if all 1 based
    int yzshift = horiz + horiz*vert;
    ebelow = num-2-yzshift;
    eabove = num+yzshift;
}

// For structured grid, find node for input to moving velocity boundary conditions. It will find the closest
// 		node to n0 but inside the object (higher number for side=-1 or lower for side=1). If
// 		position is on a node move -side nodes away.
// Also returns the gap between nodes in the input direction in dir to step through
// 		the nodes
// dir is only X_DIRECTION, Y_DIRECTION, or Z_DIRECTION, side is -1 or +1
// Handles both equal an unequal element sizes
// throws CommonException() if position has left the grid
int MeshInfo::FindShiftedNodeFromNode(int n0,double position,int dir,int side,int &nodeStep,double depth)
{
	int shift,firstNode = n0;
	bool leftGrid = false;
	
	// Find step size from inside to outside nodes in nodeStep
	// Check if left the grid
	switch(dir)
	{	case X_DIRECTION:
			nodeStep = side*xplane;
			if(position>=xmax || position<=xmin) leftGrid = true;
			break;
		case Y_DIRECTION:
			nodeStep = side*yplane;
			if(position>=ymax || position<=ymin) leftGrid = true;
			break;
		default:
			nodeStep = side*zplane;
			if(position>=zmax || position<=zmin) leftGrid = true;
			break;
	}
	
	// Exception if left the grid
	if(leftGrid)
	{   char msg[100];
		size_t msgSize=100;
		snprintf(msg,msgSize,"Moving condition starting on node %d has left the grid",n0);
		throw CommonException(msg,"");
	}

	// For equal element sizes, find first node such that
	//     depth-1 < fabs(position-firstNodes)/cell <= depth
	// In other words firstNode is last node within depth cells from the wall
	// For unequal element sizes, depth is currently 1
	// For side=1, first node is before position, for side=-1 it is after
	// For equal element algorithm, see JANOSU-14-43
	switch(dir)
	{	case X_DIRECTION:
			if(equalElementSizes)
			{	double d0 = (double)side*(position - nd[n0]->x)/grid.x;
				shift = (int)floor(depth-d0);
				if((double)shift == depth-d0-1.) shift++;
				firstNode = n0 - side*shift*xplane;
			}
			else
			{	// find firstNode before position
				if(nd[firstNode]->x<=position)
				{	firstNode += xplane;
					while(nd[firstNode]->x<=position)
						firstNode += xplane;
				}
				else
				{	while(nd[firstNode-xplane]->x>position)
						firstNode -= xplane;
				}
				// Now nd[firstNode-xplane]->x <= position < nd[firstNode]->x
				double dfirst,dcell;
				if(side==1)
				{	// For side=1 shift to get
					//    nd[firstNode]->x < position <= nd[firstNode+xplane]->x
					firstNode -= xplane;
					if(nd[firstNode]->x == position) firstNode -= xplane;
				}
				// check depth using cell width for cell more interior than one with position
				dfirst = position-nd[firstNode-side*xplane]->x;
				dcell = nd[firstNode]->x-nd[firstNode-side*xplane]->x;
				if(dfirst/dcell<depth) firstNode -= side*xplane;
			}
			break;
		
		case Y_DIRECTION:
			if(equalElementSizes)
			{	double d0 = (double)side*(position - nd[n0]->y)/grid.y;
				shift = (int)floor(depth-d0);
				if((double)shift == depth-d0-1.) shift++;
				firstNode = n0 - side*shift*yplane;
			}
			else
			{	// find firstNode before position
				if(nd[firstNode]->y<=position)
				{	firstNode += yplane;
					while(nd[firstNode]->y<=position)
						firstNode += yplane;
				}
				else
				{	while(nd[firstNode]->y>position)
						firstNode -= yplane;
				}
				// Now nd[firstNode-yplane]->y <= position < nd[firstNode]->y
				double dfirst,dcell;
				if(side==1)
				{	// For side=1 shift to get
					//    nd[firstNode]->y < position <= nd[firstNode+yplane]->y
					firstNode -= yplane;
					if(nd[firstNode]->y == position) firstNode -= yplane;
				}
				// check depth using cell width for cell more interior than one with position
				dfirst = position-nd[firstNode-side*yplane]->y;
				dcell = nd[firstNode]->y-nd[firstNode-side*yplane]->y;
				if(dfirst/dcell<depth) firstNode -= side*yplane;
			}
			break;
		
		default:
			if(equalElementSizes)
			{	double d0 = (double)side*(position - nd[n0]->z)/grid.z;
				shift = (int)floor(depth-d0);
				if((double)shift == depth-d0-1.) shift++;
				firstNode = n0 - side*shift*zplane;
			}
			else
			{	// find firstNode before position
				if(nd[firstNode]->z<=position)
				{	firstNode += zplane;
					while(nd[firstNode]->z<=position)
						firstNode += zplane;
				}
				else
				{	while(nd[firstNode]->z>position)
						firstNode -= zplane;
				}
				// Now nd[firstNode-zplane]->z <= position < nd[firstNode]->z
				double dfirst,dcell;
				if(side==1)
				{	// For side=1 shift to get
					//    nd[firstNode]->z < position <= nd[firstNode+zplane]->z
					firstNode -= zplane;
					if(nd[firstNode]->z == position) firstNode -= zplane;
				}
				// check depth using cell width for cell more interior than one with position
				dfirst = position-nd[firstNode-side*zplane]->z;
				dcell = nd[firstNode]->z-nd[firstNode-side*zplane]->z;
				if(dfirst/dcell<depth) firstNode -= side*zplane;
			}
			break;
	}
	
	return firstNode;
}

// For structured grid, last possible node for velocity gradient BC with
//		n0 as the firstNode
// dir is only X_DIRECTION, Y_DIRECTION, or Z_DIRECTION, side is -1 or +1
int MeshInfo::FindLastNode(int n0,int dir,int side)
{
	// in the rank (=0 if 2D)
	// n0 = rank*zplane + row*yplane + col + 1		// 1 based
	// where col, row, and rank are zero based
	int col = (n0-1) % yplane ;		// column for the node
	int rowStart = n0-1-col;		// node-1 at start of row containing n0 = rank*zplane + row*yplane
	int lastNode;
	if(Is3DGrid())
	{	// 3D grid
		int rank;
		switch(dir)
		{	case X_DIRECTION:
				lastNode = side<0 ? rowStart+1 : rowStart+yplane ;
				break;
			case Y_DIRECTION:
				rank = n0 % zplane;
				lastNode = side<0 ? rank*zplane+col+1 : (rank+1)*zplane-yplane+col+1 ;
			default:
				rank = n0 % zplane;
				rowStart -= rank*zplane;		// rowStart in rank 0
				lastNode = side<0 ? rowStart+col+1 : nnodes-zplane+rowStart+col+1;
				break;
		}
	}
	else
	{	// 2D gris
		switch(dir)
		{	case X_DIRECTION:
				lastNode = side<0 ? rowStart+1 : rowStart+yplane ;
				break;
			default:
				lastNode = side<0 ? col+1 : zplane-yplane+col+1 ;
				break;
		}
	}
	
	return lastNode;
}

// For structured, find element from location and return result (1-based element number)
// Calling code must be sure it is structured grid
// throws CommonException()
int MeshInfo::FindElementFromPoint(const Vector *pt,MPMBase *mptr)
{
    int theElem,col,row,zrow;
	
	if(equalElementSizes)
    {	col = (int)((pt->x-xmin)/grid.x);		// zero-based column # from 0 to horiz-1
		if(col<0 || col>=horiz)
		{	if(pt->x == xmin+horiz*grid.x)
				col = horiz-1;
			else
			{	throw CommonException("column out of range","");
			}
		}
    
		row = (int)((pt->y-ymin)/grid.y);        // zero-based row # from 0 to vert-1
		if(row<0 || row>=vert)
		{	if(pt->y == ymin+vert*grid.y)
				row = vert-1;
			else
			{	throw CommonException("row out of range","");
			}
		}
    
		// 3D
		if(grid.z > 0.)
		{   zrow = (int)((pt->z-zmin)/grid.z);   // zero-based row # from 0 to depth-1
			if(zrow<0 || zrow>=depth)
			{	if(pt->z == zmin+depth*grid.z)
					zrow = depth-1;
				else
				{	throw CommonException("rank out of range","");
				}
			}
			theElem = horiz*(zrow*vert + row) + col + 1;
		}
    
		// 2D
		else
		{   theElem = row*horiz + col + 1;
		}
	}
	
	// structured grid, but unequal element sizes
	else
    {	Vector testPt = *pt;
		
		// if has point
		if(mptr!=NULL)
		{	// check current element
			theElem = mptr->ElemID();
			if(theElements[theElem]->PtInElement(testPt))
				return theElem+1;
			
			// look at nearest neighbors
			int i=0,elemNeighbors[27];
			theElements[mptr->ElemID()]->GetListOfNeighbors(elemNeighbors);
			while(elemNeighbors[i]!=0)
			{	if(theElements[elemNeighbors[i]-1]->PtInElement(testPt))
					return elemNeighbors[i];
				i++;
			}
		}
		
		// x axis elements 0 to horiz-1
		col = BinarySearchForElement(0,pt->x,horiz-1,1);
		if(col<0)
		{	throw CommonException("column for tartan grid out of range","");
		}

		// y axis elements 0 to (vert-1)*horiz
		row = BinarySearchForElement(1,pt->y,vert-1,horiz);
		if(row<0)
		{	throw CommonException("row for tartan grid out of range","");
		}

		if(fmobj->IsThreeD())
		{	// z axis element 0 to (depth-1)*horiz*vert
			zrow = BinarySearchForElement(2,pt->z,depth-1,horiz*vert);
			if(zrow<0)
			{	throw CommonException("rank for tartan grid out of range","");
			}
			
			theElem = horiz*(zrow*vert + row) + col + 1;
		}
		
		// 2D
		else
		{   theElem = row*horiz + col + 1;
		}
	}
	
	// unstructured grid is disabled
    
    // return result
    return theElem;
}

// binary search to find element in grid with unequal element sizes
// element is >=min to <max (except <= last one's max)
int MeshInfo::BinarySearchForElement(int ax,double pt,int last,int numGap)
{
	double amin,amax;
	
	// check first element
	theElements[0]->GetRange(ax,amin,amax);
	
	// out of range
	if(pt<amin) return -1;
	
	// is it in the first element
	if(pt<amax) return 0;
	
	// check last element
	int etest = last*numGap;
	theElements[etest]->GetRange(ax,amin,amax);
	
	// out of range
	if(pt>amax) return -1;
	
	// is it in the first element
	if(pt>=amin) return last;
	
	// binary search
	int e1 = 0,e2 = last,emid;
	while(true)
	{	emid = (e1+e2)>>1;
		etest = emid*numGap;
		theElements[etest]->GetRange(ax,amin,amax);
		
		// to the left
		if(pt < amin)
			e2 = emid;
		
		// to the right
		else if(pt >= amax)
			e1 = emid;
		
		// done
		else
			break;
	}
	
	// return result, which musr be found
	return emid;
}

// For structured, find element from location and return result coordinates
// Calling code must be sure it is structured grid with equal element sizes
// col, row, and zrow are zero based
void MeshInfo::FindElementCoordinatesFromPoint(Vector *pt,int &col,int &row,int &zrow)
{
	zrow=0;
	
	if(equalElementSizes)
	{	col = (int)((pt->x-xmin)/grid.x);		// zero-based column # from 0 to horiz-1
		if(col<0 || col>=horiz)
		{	if(pt->x == xmin+horiz*grid.x)
				col = horiz-1;
			else
			{   char msg[100];
				size_t msgSize=100;
				snprintf(msg,msgSize,"column for point (%lf,%lf,%lf)",pt->x,pt->y,pt->z);
				throw CommonException(msg,"");
			}
		}
		
		row = (int)((pt->y-ymin)/grid.y);        // zero-based row # from 0 to vert-1
		if(row<0 || row>=vert)
		{	if(pt->y == ymin+vert*grid.y)
				row = vert-1;
			else
			{	char msg[100];
				size_t msgSize=100;
				snprintf(msg,msgSize,"row for point (%lf,%lf,%lf)",pt->x,pt->y,pt->z);
				throw CommonException(msg,"");
			}
		}
		
		// 3D
		if(grid.z > 0.)
		{   zrow = (int)((pt->z-zmin)/grid.z);   // zero-based row # from 0 to depth-1
			if(zrow<0 || zrow>=depth)
			{	if(pt->z == zmin+depth*grid.z)
					zrow = depth-1;
				else
				{   char msg[100];
					size_t msgSize=100;
					snprintf(msg,msgSize,"rank for point (%lf,%lf,%lf)",pt->x,pt->y,pt->z);
					throw CommonException(msg,"");
				}
			}
		}
	}
	
	else
		throw CommonException("Invalid attempt to find coordinates from point in grid with unequal element sizes","");
}

// Set custom patch grid.
// These will not be used unless xpnum*ypum*zpnum = number of processors
// zpnum is set to 1 in 2D regardless of setting here
void MeshInfo::SetCustomPatching(int xnum,int ynum,int znum)
{	xpnum = xnum;
	ypnum = ynum;
	zpnum = znum;
}

// Swt the size that is about to be set
// These will not be used unless xpnum*ypum*zpnum = number of processors
// zpnum is set to 1 in 2D regardless of setting here
void MeshInfo::SetSizeDirection(int dir)
{	sizeDir = dir;
}

// Valid and then set up one patch size direction
void MeshInfo::SetUpPatchSizes(vector<double> &psizes,int pnum)
{
	if(psizes.size()>0)
	{	// must have pnum-1 (and add last one)
		if(psizes.size()!=pnum-1)
		{	throw CommonException("Patch sizes do not match number of patches in a direction",
								  "MeshInfo::SetUpPatchSizes");
		}
		// must be ordered and less than 1
		for(int i=0;i<pnum-1;i++)
		{	if(psizes[i]<=0. || psizes[i]>=1.)
			{	throw CommonException("Patch sizes must be between 0 and 1",
									  "MeshInfo::SetUpPatchSizes");
			}
			if(i>0)
			{	if(psizes[i]<=psizes[i-1])
				{	throw CommonException("Patch sizes must montonically increase",
										  "MeshInfo::SetUpPatchSizes");
				}
			}
		}
	}
	else
	{	for(int i=1;i<pnum;i++)
			psizes.push_back(i/(double)pnum);
	}
	psizes.push_back(1.0);

}


// set from  arbitrary string data
// not thread safe due to push_back()
// throws std::bad_alloc, SAXException()
void MeshInfo::SetPatchSizes(char *bData)
{	switch(sizeDir)
	{	case 1:
			if(xsizes.size()>0)
				throw SAXException("X patch sizes set more than once.");
			if(!CommonReadHandler::GetFreeFormatNumbers(bData,xsizes,1.0))
				throw SAXException("Invalid X patch sizes block - bad number formatting");
			break;
		case 2:
			if(ysizes.size()>0)
				throw SAXException("Y patch sizes set more than once.");
			if(!CommonReadHandler::GetFreeFormatNumbers(bData,ysizes,1.0))
				throw SAXException("Invalid X patch sizes block - bad number formatting");
			break;
		default:
			if(zsizes.size()>0)
				throw SAXException("Z patch sizes set more than once.");
			if(!CommonReadHandler::GetFreeFormatNumbers(bData,zsizes,1.0))
				throw SAXException("Invalid Z patch sizes block - bad number formatting");
			break;
	}
}

// Create the patches for the grid
// Return pointer to a 0-based listed or patches[0] ... pathes[numProcs-1]
// Return NULL on memory error
// throws std::bad_alloc
GridPatch **MeshInfo::CreatePatches(int np,int numProcs)
{
	// serial or single thread is simpler
	if(numProcs<=1)
	{	return CreateOnePatch(np);
	}
    
    // try to use custom patching (if it was set)
	unsigned ndim = np==THREED_MPM ? 3 : 2 ;
	if(xpnum>0)
	{	if(ndim!=3) zpnum = 1;
		if(xpnum*ypnum*zpnum != numProcs)
		{	throw CommonException("Custom patching does not match number of processors",
								  "MeshInfo::CreatePatches");
		}
	}
	
	// If not custom, arrange by prime factors
	// If unequal numbers, put more in longer directions
    if(xpnum<0)
	{	// get prime factors in ascending order
		vector<int> factors;
		PrimeFactors(numProcs,factors);
		while(factors.size()<ndim) factors.push_back(1);
		while(factors.size()>ndim)
		{	std::sort(factors.begin(),factors.end());
			int newFactor = factors[0]*factors[1];
			factors.erase(factors.begin());
			factors[0]=newFactor;
		}
		std::sort(factors.begin(),factors.end());
    
		// convert to numbers of elements in each direction and use larger
        // prime factors in the longer directions
        if(ndim==2)
        {	if(horiz>=vert)
            {	// x >= y
                xpnum = factors[1];
                ypnum = factors[0];
            }
            else
            {	// y > x
                ypnum = factors[1];
                xpnum = factors[0];
            }
            zpnum=1;
        }
        else
        {   if(horiz>=vert && horiz>=depth)
            {	//	x >= y and z
                xpnum = factors[2];
                if(vert>=depth)
                {	// x >= y >= z
                    ypnum = factors[1];
                    zpnum = factors[0];
                }
                else
                {	// x >= z > y
                    zpnum = factors[1];
                    ypnum = factors[0];
                }
            }
            else if(vert>=horiz && vert>=depth)
            {	// y >= x and z
                ypnum = factors[2];
                if(horiz>=depth)
                {	// y >= x >= z
                    xpnum = factors[1];
                    zpnum = factors[0];
                }
                else
                {	// y >= z > x
                    zpnum = factors[1];
                    xpnum = factors[0];
                }
            }
            else
            {	// z > x and z
                zpnum = factors[2];
                if(horiz>=vert)
                {	// z > x >= y
                    xpnum = factors[1];
                    ypnum = factors[0];
                }
                else
                {	// z > y > x
                    ypnum = factors[1];
                    xpnum = factors[0];
                }
            }
        }
    }
	
    // get patch sizes (used only when regularly spaced)
	xPatchSize = max(int(horiz/xpnum+.5),1);
	yPatchSize = max(int(vert/ypnum+.5),1);
    if(ndim==2)
        zPatchSize = 1;
    else
        zPatchSize = max(int(depth/zpnum+.5),1);
    
    // alloc space for patches - exit on memory error
    int totalPatches = xpnum*ypnum*zpnum;
	GridPatch **patch = new (nothrow) GridPatch *[totalPatches];
    if(patch==NULL) return NULL;

	// create the patches
	int pnum=0;
	
	// optionally use tartan-style patch grid (if any direction was set)
	if(xsizes.size()>0 || ysizes.size()>0 || (ndim==3 && zsizes.size()>0))
	{	// tartan patches are requested
		tartanPatches = true;
		
		// validate and/or create sizes
		SetUpPatchSizes(xsizes,xpnum);
		SetUpPatchSizes(ysizes,ypnum);
		if(ndim==3)
			SetUpPatchSizes(zsizes,zpnum);
		else
		{	zsizes.clear();
			zsizes.push_back(1.0);
		}
	}
	
	// Create tartan patches if requested
	if(tartanPatches)
	{	// allocate pointers
		xtp = (int *)malloc(sizeof(int)*xpnum);
		ytp = (int *)malloc(sizeof(int)*ypnum);
		ztp = (int *)malloc(sizeof(int)*zpnum);
		
		// set break points and create patches
		int x1,x2,y1,y2,z1=1,z2;
		for(int k=0;k<zpnum;k++)
		{	z2 = int(depth*zsizes[k]+0.001);		// for 2D z1=1 and z2=depth=0
			ztp[k] = z2;
			if(z2<=z1 && ndim==3)
			{	throw CommonException("A Z patch block has zero elements",
									  "MeshInfo::CreatePatches");
			}
			y1 = 1;
			for(int j=0;j<ypnum;j++)
			{	y2 = int(vert*ysizes[j]+.001);
				ytp[j] = y2;
				if(y2<=y1)
				{	throw CommonException("A Y patch block has zero elements",
										  "MeshInfo::CreatePatches");
				}
				x1 = 1;
				for(int i=0;i<xpnum;i++)
				{	x2 = int(horiz*xsizes[i]+.001);
					xtp[i] = x2;
					if(x2<=x1)
					{	throw CommonException("An X patch block has zero elements",
											  "MeshInfo::CreatePatches");
					}

					// patch x1 to x2 and y1 to y2 (1 based)
					patch[pnum] = new GridPatch(x1,x2,y1,y2,z1,z2);
					if(!patch[pnum]->CreateGhostNodes())
					{	delete [] patch;
						return NULL;
					}
					pnum++;
					
					x1 = x2+1;
				}
				y1 = y2+1;
			}
			z1 = z2+1;
		}
	}
	
	// If tartan patch grid not set, create reqgular patch grid now
	if(pnum==0)
	{	int x1,x2,y1,y2,z1=1,z2;
		for(int k=1;k<=zpnum;k++)
		{	z2 = k==zpnum ? depth : z1+zPatchSize-1;		// for 2D z1=1 and z2=depth=0
			y1 = 1;
			for(int j=1;j<=ypnum;j++)
			{	y2 = j==ypnum ? vert : y1+yPatchSize-1;
				x1 = 1;
				for(int i=1;i<=xpnum;i++)
				{	x2 = i==xpnum ? horiz : x1+xPatchSize-1;
					
					// patch x1 to x2 and y1 to y2 (1 based)
					patch[pnum] = new GridPatch(x1,x2,y1,y2,z1,z2);
					if(!patch[pnum]->CreateGhostNodes())
					{	delete [] patch;
						return NULL;
					}
					pnum++;
					
					x1 = x2+1;
				}
				y1 = y2+1;
			}
			z1 = z2+1;
		}
	}
	
	// fill patches with particles, but put reservoir particles in the reservoir
	int pn;
	for(int p=0;p<nmpms;p++)
	{	if(mpm[p]->InReservoir() && mpmReservoir!=NULL)
			mpmReservoir->AddParticle(mpm[p]);
	   	else
		{	pn = GetPatchForElement(mpm[p]->ElemID());
			if(pn<0 || pn>=totalPatches)
			{	delete [] patch;
				return NULL;
			}
			patch[pn]->AddParticle(mpm[p]);
		}
	}

    // return array of patches
    return patch;
}

// Create a single patch for the grid and patch has no ghost nodes
// Return pointer to a 0-based listed or patches[0]
// Return NULL on memory error
// throws std::bad_alloc
GridPatch **MeshInfo::CreateOnePatch(int np)
{
	// a single patch
	xpnum = ypnum = zpnum = 1;
	
	xPatchSize = horiz;
	yPatchSize = vert;
	zPatchSize = np==THREED_MPM ? depth : 1 ;
    
    // alloc space for patches - exit on memory error
	GridPatch **patch = new (nothrow) GridPatch *[1];
    if(patch==NULL) return NULL;
	
	// one patch but no ghost nodes
	patch[0] = new GridPatch(1,xPatchSize,1,yPatchSize,1,np==THREED_MPM ? depth : 0);
	
	// fill patch with all particles, except put reservoir particles in the reservoir
	for(int p=0;p<nmpms;p++)
	{	if(mpm[p]->InReservoir() && mpmReservoir!=NULL)
			mpmReservoir->AddParticle(mpm[p]);
		else
			patch[0]->AddParticle(mpm[p]);
	}
    
    // return array with the one patch
    return patch;
}

#pragma mark MeshInfo::Multimaterial contact info

// Called if multimaterial mode active: Print mode settings and contact law details
// throws CommonException()
void MeshInfo::MaterialOutput(void)
{
	// Global material contact law (must be set, if not force to frictionless)
	materialContactLawID = MaterialBase::GetContactLawNum(materialContactLawID);
	if(materialContactLawID<0)
		throw CommonException("Multimaterial mode must select a default contact law","MeshInfo::MaterialOutput");
	materialContactLaw = (ContactLaw *)theMaterials[materialContactLawID];
	cout << "Default Contact Law: " << materialContactLaw->name << " (number " << (materialContactLawID+1) << ")" << endl;
	if(materialContactLaw->IsImperfectInterface()) hasImperfectInterface=true;
	
	// request volume gradient extrapolations
	volumeGradientIndex = 0;
	
	// print contact detection method
	cout << "Contact Detection: (Normal dv < 0) & (Normal cod < 0)" << endl;
	OutputContactByDisplacements(materialNormalMethod>=LINEAR_REGRESSION,contactByDisplacements,positionCutoff);
	cout << "Normal Calculation: ";
	switch(materialNormalMethod)
	{	case MAXIMUM_VOLUME_GRADIENT:
			cout <<                     "Gradient of material or paired material (if all nonrigid), or the rigid" << endl;
			cout << "                    material (if one rigid material), that has highest magnitude. When has" << endl;
			cout << "                    rigid material, prefer rigid material gradient with bias factor = " << rigidGradientBias;
			break;
		case MAXIMUM_VOLUME:
			cout << " gradient of material with maximum volume";
			break;
		case AVERAGE_MAT_VOLUME_GRADIENTS:
			cout <<                     "Volume-weighted mean gradient of material and other materials lumped (if all" << endl;
			cout << "                    nonrigid), on just the rigid material (if one rigid material). When has" << endl;
			cout << "                    rigid material, prefer rigid material gradient with bias factor = " << rigidGradientBias;
			break;
		case EACH_MATERIALS_MASS_GRADIENT:
			cout << "Each material's own mass gradient";
			break;
		case SPECIFIED_NORMAL:
			volumeGradientIndex = -1;
			cout << "Use the specified normal of ";
			PrintVector("",&contactNormal);
			break;
		case LINEAR_REGRESSION:
			cout <<                     "Use linear regression";
			//if(rigidGradientBias>=10. && (nmpmsRC>nmpmsNR))
			if(rigidGradientBias>=10.)
				cout << " but switch to rigid material gradient for contact with rigid materials";
			else
			{	rigidGradientBias = 1.;
				volumeGradientIndex = -1;
			}
			cout << ".";
			break;
		case LOGISTIC_REGRESSION:
			cout <<                     "Use logistic regression";
			//if(rigidGradientBias>=10. && (nmpmsRC>nmpmsNR))
			if(rigidGradientBias>=10.)
				cout << " but switch to rigid material gradient for contact with rigid materials";
			else
			{	rigidGradientBias = 1.;
				volumeGradientIndex = -1;
			}
			cout << ".";
			warnLRConvergence = warnings.CreateWarning("Logistic regression did not converge and found negative separation",-1,0);
			break;
		default:
			break;
	}
	rigidGradientBias*=rigidGradientBias;       // algorithms assume it is squared (to compare to squared mag of other vectors)
	cout << endl;
	
	// lumping method
	cout << "3+ Material Contact Nodes: ";
	switch(lumpingMethod)
	{
#ifdef THREE_MAT_CONTACT
		case EXPLICIT_PAIRS:
			cout << "Explicitly handle one or two pairs in contact";
			break;
#endif // end THREE_MAT_CONTACT
		case LUMP_OTHER_MATERIALS:
		default:
			lumpingMethod = LUMP_OTHER_MATERIALS;
			cout << "Lump other materials into a virtual material";
			break;
	}
	cout << endl;
}

// return true if normals methods using regression methods and needs
// This alerts code to build list of material points seeing contact node
bool MeshInfo::UsingRegressionForContact(void)
{	if((materialNormalMethod == LINEAR_REGRESSION) ||
	   (materialNormalMethod == LOGISTIC_REGRESSION))
	{	return true;
	}
	return false;
}

// set contact normal when normal is specified
void MeshInfo::SetContactNormal(double polarAngle,double aximuthAngle)
{
	double angle,sinp;
	
	if(fmobj->IsThreeD())
	{	angle = PI_CONSTANT*polarAngle/180.;
		contactNormal.z = cos(angle);
		sinp = sin(angle);
	}
	else
	{	contactNormal.z = 0.;
		sinp = 1.0;
	}
	angle = PI_CONSTANT*aximuthAngle/180.;
	contactNormal.x = cos(angle)*sinp;
	contactNormal.y = sin(angle)*sinp;
}

// prepare array for material contact details
// throws CommonException()
void MeshInfo::MaterialContactPairs(int maxFields)
{
	// create double array of pairs
	mmContactLaw = new (nothrow) ContactLaw **[maxFields];
	if(mmContactLaw==NULL)
	{	throw CommonException("Memory error creating contact pairs array","MeshInfo::MaterialContactPairs");
	}
	
	// fill all pairs with default material properties
	int i,j;
	for(i=0;i<maxFields-1;i++)
	{	mmContactLaw[i] = new (nothrow) ContactLaw *[maxFields-1-i];
		if(mmContactLaw[i]==NULL)
		{	throw CommonException("Memory error creating contact pairs array","MeshInfo::MaterialContactPairs");
		}
		
		// to default law
		for(j=i+1;j<maxFields;j++) mmContactLaw[i][j-i-1] = mpmgrid.materialContactLaw;
	}
	
	// check all active materials and change laws that were specified
	for(i=0;i<nmat;i++)
	{	int mati=theMaterials[i]->GetField();			// may be a shared field
		if(mati<0) continue;							// skip if not used
		
		// loop over all other materials
		for(j=0;j<nmat;j++)
		{	int matj=theMaterials[j]->GetField();		// may be a shared field
			if(matj<0 || i==j) continue;				// skip if no field or same material
			
			// look from custom friction from mat i to mat j
			int pairContactID=theMaterials[i]->GetContactToMaterial(j+1);
			if(pairContactID<0) continue;
			pairContactID = MaterialBase::GetContactLawNum(pairContactID);
			
			// setting more than one shared material overwrite previous ones
			if(mati<matj)
			{	mmContactLaw[mati][matj-mati-1]=(ContactLaw *)theMaterials[pairContactID];
				if(mmContactLaw[mati][matj-mati-1]->IsImperfectInterface()) mpmgrid.hasImperfectInterface=true;
			}
			else
			{	mmContactLaw[matj][mati-matj-1]=(ContactLaw *)theMaterials[pairContactID];
				if(mmContactLaw[matj][mati-matj-1]->IsImperfectInterface()) mpmgrid.hasImperfectInterface=true;
			}
		}
	}
}

// material contact law for field mati to field matj
ContactLaw *MeshInfo::GetMaterialContactLaw(int mati,int matj)
{	// index based on smaller of the two indices
	return mati<matj ? mmContactLaw[mati][matj-mati-1] : mmContactLaw[matj][mati-matj-1] ;
}

#pragma mark MeshInfo:Accessors

// find data for creating symmetry BCs
// input axis (x,y,z), sym location, and direction (-1 or +1)
// output node number (from 1 to count), element size outside and inside the node (gridout,gridin)
bool MeshInfo::FindMeshLineForBCs(int axis,double gridsym,int symdir,int &nnum,double &gridout,double &gridin)
{
	int gap,count;
	
	// for each axis, make settings and check if outside the grid
	switch(axis)
	{	case X_DIRECTION:
			gap = xplane;
			count = horiz+1;
			if((symdir<0 && gridsym<xmin) || (symdir>0 && gridsym>xmax))
				return false;
			break;
			
		case Y_DIRECTION:
			gap = yplane;
			count = vert+1;
			if((symdir<0 && gridsym<ymin) || (symdir>0 && gridsym>ymax))
				return false;
			break;
			
		default:
			gap = zplane;
			count = depth+1;
			if((symdir<0 && gridsym<zmin) || (symdir>0 && gridsym>zmax))
				return false;
			break;
	}
	
	// switch if from maximum end
	int firstNode=1;
	int lastNode = firstNode + (count-1)*gap;
	if(symdir>0)
	{	int tmp = firstNode;
		firstNode = lastNode;
		lastNode = tmp;
		gap = -gap;
	}
	
	// step to the node
	nnum = 1;
	int j = firstNode;
	double pt,prevpt=0.,nextpt=0.;
	while(nnum <= count)
	{	if(axis==X_DIRECTION)
		{	pt = nd[j]->x;
			if(j!=lastNode) nextpt = nd[j+gap]->x;
		}
		else if(axis==Y_DIRECTION)
		{	pt = nd[j]->y;
			if(j!=lastNode) nextpt = nd[j+gap]->y;
		}
		else
		{	pt = nd[j]->z;
			if(j!=lastNode) nextpt = nd[j+gap]->z;
		}
		
		// check if close enough (closer than grid tolerance, which is 10% of cell)
		if(fabs(pt-gridsym)<ElementBase::gridTolerance/100.)
		{	// fails if node is on an edge
			if(j==firstNode || j==lastNode) return false;
			gridout = fabs(pt-prevpt);
			gridin = fabs(nextpt-pt);
			if(symdir>0) nnum = count+1-nnum;
			break;
		}
		
		// if has passed the desired position, then not on a grid line
		if((symdir<0 && gridsym<pt) || (symdir>0 && gridsym>pt))
			return false;
		
		// next node
		j += gap;
		prevpt = pt;
		nnum++;
	}
	
	// break above is correct result
	return nnum<=count ? true : false ;
}

// given zero based element number, return zero-based patch number
int MeshInfo::GetPatchForElement(int iel)
{
	int row,col,rank;
	int prow,pcol,prank;
	
	if(depth>0)
	{	// 3D
		int perSlice = horiz*vert;		// number in each slice
		rank = iel/perSlice;			// zero based
		int snum = iel % perSlice;		// number in slice (0 to horiz*vert-1)
		col = snum % horiz;				// col 0 to horiz-1
		row = snum/horiz;				// zero based
		if(tartanPatches)
		{	pcol = 0;
			while(col>=xtp[pcol]) pcol++;
			prow = 0;
			while(row>=ytp[prow]) prow++;
			prank = 0;
			while(rank>=ztp[prank]) prank++;
		}
		else
		{	pcol = min(col/xPatchSize,xpnum-1);
			prow = min(row/yPatchSize,ypnum-1);
			prank = min(rank/zPatchSize,zpnum-1);
		}
		return xpnum*(ypnum*prank + prow) + pcol;
	}
	
	// 2D
	col = iel % horiz;			// col 0 to horiz-1
	row = iel/horiz;			// zero based
	if(tartanPatches)
	{	pcol = 0;
		while(col>=xtp[pcol]) pcol++;
		prow = 0;
		while(row>=ytp[prow]) prow++;
	}
	else
	{	pcol = min(col/xPatchSize,xpnum-1);
		prow = min(row/yPatchSize,ypnum-1);
	}
	return xpnum*prow + pcol;
}

// given zero based element number, find row, col, rank
// row, col, and rank are zero based
// Warning: assume structure grid (and not checked here)
void MeshInfo::GetElementsRCR(int iel,int &row,int &col,int &rank)
{
    if(depth>0)
    {    // 3D
        int perSlice = horiz*vert;        // number in each slice
        rank = iel/perSlice;            // zero based
        int snum = iel % perSlice;        // number in slice (0 to horiz*vert-1)
        col = snum % horiz;                // col 0 to horiz-1
        row = snum/horiz;                // zero based
        return;
    }
    
    // 2D
    col = iel % horiz;            // col 0 to horiz-1
    row = iel/horiz;            // zero based
    return;
}

// set grid style (zcell=0 if 2D grid)
// Options:
//	  style=NOT_CARTESIAN means not aligned with x,y,z axes
//		just set cartesian to NOT_CARTESIAN (or _3D) and equalElementSizes to false
//    style=VARIABLE_ORTHOGONAL_GRID
//	  style=VARIABLE_RECTANGULAR_GRID
//      3D or 2D grid with unequal element sizes, set equalElementSizes to false
//		grid.x=grid.y=0 and grid.z=0
//    style=SQUARE_GRID on input means all elements the same size as subsets, but adjust as follows
//      SQUARE_GRID = 2D and dx = dy
//      CUBIC_GRID = 3D and dx = dy = dz
//      RECTANGULAR_GRID = 2D and dx != dy
//      ORTHOGONAL_GRID = 3D and one of dx, dy, dz != others
//		set equalElementSizes to true
void MeshInfo::SetCartesian(int style,double xcell,double ycell,double zcell)
{
	grid.x=xcell;
	grid.y=ycell;
	grid.z=zcell;
	cartesian=style;
	equalElementSizes = false;
	
	// if variable, then done
    if(style == VARIABLE_ORTHOGONAL_GRID || style == VARIABLE_RECTANGULAR_GRID)
    {   return;
    }
	
	// if not cartensian, set 2D or 2D
	else if(style==NOT_CARTESIAN)
	{	if(fmobj->IsThreeD())
			cartesian = NOT_CARTESIAN_3D;
	}
    
	else if(style>0)
	{	grid.x=xcell;
		grid.y=ycell;
		grid.z=zcell;
		if(DbleEqual(xcell,ycell))
		{	if(DbleEqual(zcell,0.))
                cartesian=SQUARE_GRID;      // 2D, xcell=ycell
            else if(DbleEqual(xcell,zcell))
                cartesian=CUBIC_GRID;       // 3D, xcell=ycell=zcell
            else
                cartesian=ORTHOGONAL_GRID;  // 3D, xcell=ycell != zcell
		}
		else if(DbleEqual(zcell,0.))
			cartesian=RECTANGULAR_GRID;     // 2D, xcell != ycell
		else
			cartesian=ORTHOGONAL_GRID;      // 3D, xcell != ycell, zcell not checked
		
		// equal element sizes
		equalElementSizes = true;
	}
}

// set grid style (only set by <Grid> command). For 2D, d=0 and z is the specified thickness, or 1.0 by default,
//   or 1.0 if axisymmetric (so cell Volume is area)
// horiz is number of elements in x direction, vert is y direction, depth in z direction (or 0 in 2D)
//		cellVolume found and avgCellSize is average length of cell edges
// Also called for structured grid with variable elements
//		cellVolume = avgCellSize = 0 (because grid was set to zero)
// Always called by any structured grid, which current is all problems
void MeshInfo::SetElements(int h,int v,int d,double x,double y,double z,double xmx,double ymx,double zmx)
{	horiz=h;
	vert=v;
	depth=d;
	xmin=x;
	ymin=y;
	zmin=z;
	xmax=xmx;
	ymax=ymx;
	zmax=zmx;
	
	// spacing between nodes numbers along each axis
	xplane=1;
	yplane=horiz+1;
	zplane=(horiz+1)*(vert+1);
	if(depth>0)
	{	// 3D grid
		totalElems=horiz*vert*depth;
		cellVolume=grid.x*grid.y*grid.z;
        avgCellSize=(grid.x+grid.y+grid.z)/3.;
	}
	else
	{	// 2D grid
		totalElems=horiz*vert;
		cellVolume=grid.x*grid.y*zmin;
        avgCellSize=(grid.x+grid.y)/2.;
	}
}

// return cartesian setting
int MeshInfo::GetCartesian(void) { return cartesian; }

// see if 3D grid
bool MeshInfo::Is3DGrid(void) { return cartesian > BEGIN_3D_GRIDS; }

// find element number in lower left hand corner of a structure
// grid, but not one on the edge
int MeshInfo::GetCornerNonEdgeElementNumber(void)
{   return Is3DGrid() ? 2+horiz+horiz*vert : 2+horiz ;
}

// check if mesh allowed for axisymmetric GIMP shape functions
// Need equal element sizes and if min<2*cell, need 0,cell,... on grid lines
bool MeshInfo::ValidateASForGIMP(void)
{	if(!equalElementSizes) return false;
	if(xmin>=2.*grid.x) return true;
	
	// verify r=cell on grid line
	double ccell = 1.-xmin/grid.x;
	if(!DbleEqual(ccell,(double)int(ccell)) && !DbleEqual(ccell,(double)(int(ccell)+1)))
	    return false;
	return true;
}

// return minimum cell dimension if the grid is cartesian (-1 if not)
// should find the minimum size in preliminary calculations (after setting cartesian)
//     and store result for later use
double MeshInfo::GetMinCellDimension(void)
{
	if(cellMinSize>0.) return cellMinSize;
	
	switch(cartesian)
	{	case SQUARE_GRID:
		case CUBIC_GRID:
			cellMinSize=grid.x;
			break;
			
		case RECTANGULAR_GRID:
			cellMinSize=fmin(grid.x,grid.y);
			break;
			
		case ORTHOGONAL_GRID:
			cellMinSize=fmin(grid.x,grid.y);
			cellMinSize=fmin(cellMinSize,grid.z);
			break;
			
		default:
			break;
	}
	
	// search all elements
	if(cellMinSize<0.)
	{	int i;
		ElementBase *elem;
		cellMinSize=1.e30;
		for(i=0;i<nelems;i++)
		{	elem=theElements[i];
			cellMinSize = fmin(cellMinSize,elem->GetDeltaX());
			cellMinSize = fmin(cellMinSize,elem->GetDeltaY());
			if(fmobj->IsThreeD())
				cellMinSize = fmin(cellMinSize,elem->GetDeltaZ());
		}
	}
	
	// return result that is now saved
	return cellMinSize;
}

// when create particle, keep track of the smallest particle
void MeshInfo::TrackMinParticleSize(Vector newParticleSize)
{	minParticleSize.x = fmin(newParticleSize.x,minParticleSize.x);
	minParticleSize.y = fmin(newParticleSize.y,minParticleSize.y);
	minParticleSize.z = fmin(newParticleSize.z,minParticleSize.z);
}

// return vector for smallest particle size (z is thickness in 2D)
Vector MeshInfo::GetGlobalMinParticleSize(void) const { return minParticleSize; }

// return shortest particle length
double MeshInfo::GetGlobalMinParticleLength(void) const
{	double minLength = fmin(minParticleSize.x,minParticleSize.y);
	if(depth>0)
	{	// check z for 3D
		minLength = fmin(minLength,minParticleSize.z);
	}
	return 2.*minLength;
}

// horiz will be set >0 if generated by <Grid> command (but elements may have different sizes)
// currently only structured grids are allowed
bool MeshInfo::IsStructuredGrid(void) const { return horiz>0; }

// horiz will be set >0 if generated by <Grid> command and won't have variable element sizes
bool MeshInfo::IsStructuredEqualElementsGrid(void) const { return equalElementSizes; }

// number of points in the three directions (1 more than number of elements, ptz will be 1 if 2D)
// Feature that calls this method must require the problem to have a structured <Grid>
void MeshInfo::GetGridPoints(int *ptx,int *pty,int *ptz)
{	*ptx=horiz+1;
	*pty=vert+1;
	*ptz=depth+1;
}

// cell volume in mm^3
// Feature that calls this method must require the problem to have a structured <Grid>
double MeshInfo::GetCellVolume(NodalPoint *ndptr)
{
	// structured grid with equal element sizes
	if(equalElementSizes) return cellVolume;
	
	// structured grid (assumed) with variable element sizes
	double dx1,dx2,dy1,dy2,dz1,dz2;
	GetLocalCellSizes(ndptr,dx1,dx2,dy1,dy2,dz1,dz2);
	return 0.125*(dx1+dx2)*(dy1+dy2)*(dz1+dz2) ;
}

// get element lengths around 1 node in variable element cell.
// if on edge, one direction will be zero.
void MeshInfo::GetLocalCellSizes(NodalPoint *ndptr,double &dx1,double &dx2,double &dy1,double &dy2,double &dz1,double &dz2)
{	// 1 base node number, which is how they are stored in nd[]
	int num = ndptr->num;
	
	// x axis
	double x = ndptr->x;
	dx1 = num<2 ? 0. : x - nd[num-1]->x;
	if(dx1<0.) dx1 = 0.;
	dx2 = num>=nnodes ? 0. : nd[num+1]->x - x;
	if(dx2<0.) dx2 = 0.;
	
	// y axis
	double y = ndptr->y;
	dy1 = num-yplane<1 ? 0 : y - nd[num-yplane]->y;
	if(dy1<0.) dy1 = 0.;
	dy2 = num+yplane>nnodes ? 0 : nd[num+yplane]->y - y;
	if(dy2<0.) dy2 = 0.;
	
	// z axis
	if(depth>0)
	{	double z = ndptr->z;
		dz1 = num-zplane<1 ? 0 : z - nd[num-zplane]->z;
		if(dz1<0.) dz1 = 0.;
		dz2 = num+zplane>nnodes ? 0 : nd[num+zplane]->z - z;
		if(dz2<0.) dz2 = 0.;
	}
	else
	{	dz1 = zmin;
		dz2 = zmin;
	}
}

// Get ratio of element size to the right to element size to the left
//		if mirrorSpacing>0 get left to right ratio instead
// Used by rigid particle mirrored BCs and must be an interior node
double MeshInfo::GetCellRatio(NodalPoint *ndptr,int dir,int mirrorSpacing)
{
	// easy if regular grid
	if(equalElementSizes) return 1.;
	
	// rest for Tartan grid
	int num = ndptr->num;
	double dx1,dx2;
	if(dir==X_DIRECTION)
	{	// x axis
		double x = ndptr->x;
		dx1 = x - nd[num-1]->x;
		dx2 = nd[num+1]->x - x;
	}
	
	else if(dir==Y_DIRECTION)
	{	// y axis
		double y = ndptr->y;
		dx1 = y - nd[num-yplane]->y;
		dx2 = nd[num+yplane]->y - y;
	}
	
	else
	{	// z axis
		double z = ndptr->z;
		dx1 = z - nd[num-zplane]->z;
		dx2 = nd[num+zplane]->z - z;
	}
	
	return mirrorSpacing>0 ? dx1/dx2 : dx2/dx1 ;
}

// cell size in mm
// Feature that calls this method must require the problem to have a structured <Grid>
double MeshInfo::GetAverageCellSize(MPMBase *mptr)
{
	// precalculated for equal element sizes
	if(equalElementSizes) return avgCellSize;
	
	// any other grid get from the element
	Vector box = theElements[mptr->ElemID()]->GetDeltaBox();
	if(Is3DGrid())
		return (box.x+box.y+box.z)/3.;
	else
		return 0.5*(box.x+box.y);
}

// grid thickness. For 3D returns z extent but not used in 3D
// assume structured grid
double MeshInfo::GetThickness(void)
{	return depth>0 ? zmax-zmin : zmin;
}

// get thickness or 1.0 if not a structured grid
double MeshInfo::GetDefaultThickness()
{	double gthick=GetThickness();
	return gthick>0. ? gthick : 1.0 ;
}

// find hperp distance used in contact calculations in interface force calculations
Vector MeshInfo::GetPerpendicularDistance(Vector *norm,NodalPoint *ndptr)
{
	// magnitude of hperp and hperp before and after the node
	// correct for equal elements and cubic or square grids
	Vector dist = MakeVector(grid.x, 1., 1.);
    
    // Angled path correction method 1: hperp  is distance to ellipsoid through cell corners
    //    defined by tangent vector. In 3D, also multiply by distance to ellipsoid along
    //    n X t (which is along z axis for 2D)
    // In 2D and 3D the dist is equal to grid spacing if grid.x=grid.y=grid.z and therefore this
    //    whole block gets skipped
    // See JANOSU-6-60 and JANOSU-6-74 and method #1 in paper
    if(Is3DGrid())
	{
		if(!IsStructuredEqualElementsGrid())
		{	// Tartan 3D Grid
			
			// get element sizes
			double dx1,dx2,dy1,dy2,dz1,dz2;
			GetLocalCellSizes(ndptr,dx1,dx2,dy1,dy2,dz1,dz2);
			
			// get two tangents
			Vector t1,t2;
			if(norm->z>norm->x && norm->z>norm->y)
			{	// z is largest
				t1 = MakeVector(0.,-norm->z,norm->y);
			}
			else
			{	// x  or yis largest)
				t1 = MakeVector(-norm->y,norm->x,0.);
			}
			
			// second tangent from t1 X norm
			t2.x = norm->y*t1.z - norm->z*t1.y;
			t2.y = norm->z*t1.x - norm->x*t1.z;
			t2.z = norm->x*t1.y - norm->y*t1.x;
			
			// first tangent
			double a1,a2,b1,b2,c1,c2;
			ScaleTangent(t1.x,dx1,dx2,a1,a2);
			ScaleTangent(t1.y,dy1,dy2,b1,b2);
			ScaleTangent(t1.z,dz1,dz2,c1,c2);
			double h11 = 1./sqrt(a1*a1 + b1*b1 +c1*c1);
			double h12 = 1./sqrt(a2*a2 + b2*b2 + c2*c2);
			
			ScaleTangent(t2.x,dx1,dx2,a1,a2);
			ScaleTangent(t2.y,dy1,dy2,b1,b2);
			ScaleTangent(t2.z,dz1,dz2,c1,c2);
			double h21 = 1./sqrt(a1*a1 + b1*b1 + c1*c1);
			double h22 = 1./sqrt(a2*a2 + b2*b2 + c2*c2);
			
			dist.x = 0.5*(dx1+dx2)*(dy1+dy2)*(dz1+dz2)/((h11+h12)*(h21+h22));

			// partition two sides of the node
			ScaleTangent(norm->x,dx1,dx2,a1,a2);
			ScaleTangent(norm->y,dy1,dy2,b1,b2);
			ScaleTangent(norm->z,dz1,dz2,c1,c2);
			double rn = sqrt((a1*a1 + b1*b1 + c1*c1)/(a2*a2 + b2*b2 + c2*c2));
			dist.y = 0.5*(1.+rn);				// hperp/h1
			dist.z = dist.y/rn;					// hperp/h2
		}
        else if(cartesian!=CUBIC_GRID)
		{	// get two tangents
			Vector t1,t2;
			if(norm->z>norm->x && norm->z>norm->y)
			{	// z is largest
				t1 = MakeVector(0.,-norm->z,norm->y);
			}
			else
			{	// x  or yis largest)
				t1 = MakeVector(-norm->y,norm->x,0.);
			}
			
			// second tangent from t1 X norm
			t2.x = norm->y*t1.z - norm->z*t1.y;
			t2.y = norm->z*t1.x - norm->x*t1.z;
			t2.z = norm->x*t1.y - norm->y*t1.x;
			
			double a1 = t1.x/grid.x;
			double b1 = t1.y/grid.y;
			double c1 = t1.z/grid.z;
			double a2 = t2.x/grid.x;
			double b2 = t2.y/grid.y;
			double c2 = t2.z/grid.z;
			dist.x = grid.x*grid.y*grid.z*sqrt((a1*a1 + b1*b1 + c1*c1)*(a2*a2 + b2*b2 + c2*c2));
        }
    }
	else if(!IsStructuredEqualElementsGrid())
	{	// Tartan 2D grid
		
		// get element sizes
		double dx1,dx2,dy1,dy2,dz1,dz2;
		GetLocalCellSizes(ndptr,dx1,dx2,dy1,dy2,dz1,dz2);
		
		// if locally regular
		if(DbleEqual(dx1,dx2) && DbleEqual(dy1,dy2))
		{   double a=0.5*(dx1+dx2)*norm->x;
			double b=0.5*(dy1+dy2)*norm->y;
			dist.x = sqrt(a*a + b*b);
		}
		
		else
		{	// tangent Vector
			Vector t1 = MakeVector(-norm->y, norm->x, 0.);
		
			// get size to ellipse in tangential direction
			double a1,a2,b1,b2;
			ScaleTangent(t1.x,dx1,dx2,a1,a2);
			ScaleTangent(t1.y,dy1,dy2,b1,b2);
			double h1 = 1./sqrt(a1*a1 + b1*b1);
			double h2 = 1./sqrt(a2*a2 + b2*b2);
			dist.x = 0.5*(dx1+dx2)*(dy1+dy2)/(h1+h2);
			
			// partition two sides of the node
			ScaleTangent(norm->x,dx1,dx2,a1,a2);
			ScaleTangent(norm->y,dy1,dy2,b1,b2);
			double rn = sqrt((a1*a1 + b1*b1)/(a2*a2 + b2*b2));
			dist.y = 0.5*(1.+rn);				// hperp/h1
			dist.z = dist.y/rn;					// hperp/h2
		}
	}
    else if(cartesian!=SQUARE_GRID)
    {   double a=grid.x*norm->x;
        double b=grid.y*norm->y;
        dist.x = sqrt(a*a + b*b);
    }

    return dist;
}
	
// scale tangent vector
// a1 is term in opposite direction as tc and a2 in tc direction
void MeshInfo::ScaleTangent(double tc,double d1,double d2,double &a1,double &a2)
{	if(tc>0.)
	{	a1 = tc/d1;
		a2 = tc/d2;
	}
	else
	{	a1 = tc/d2;
		a2 = tc/d1;
	}
}

// Cell size for grid with constant element size only
Vector MeshInfo::GetCellSize(void) { return grid; }
