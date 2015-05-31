/*********************************************************************
    MeshInfo.cpp
    Nairn Research Group MPM Code
    
    Created by John Nairn on 5/16/06.
    Copyright (c) 2006, All rights reserved.
	
	Dependencies
		none
*********************************************************************/

#include "NairnMPM_Class/MeshInfo.hpp"
#include "Patches/GridPatch.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Elements/ElementBase.hpp"
#include <algorithm>

// global class for grid information
MeshInfo mpmgrid;

#pragma mark MeshInfo:Constructors and Destructors

MeshInfo::MeshInfo(void)
{
	cartesian=UNKNOWN_GRID;
	horiz=0;                        // also flag that used <Grid> command (i.e. structured grid)
	contactByDisplacements=TRUE;	// contact by displacements
	positionCutoff=0.8;             // element fraction when contact by positions
	cellMinSize=-1.;				// calculated when needed
}

#pragma mark MeshInfo:Methods

// output grid info
void MeshInfo::Output(int pointsPerCell,bool isAxisym)
{
    char fline[200];

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
                sprintf(fline," orthogonal grid with dR: %g dZ: %g",grid.x,grid.y);
            else
                sprintf(fline," orthogonal grid with dx: %g dy: %g",grid.x,grid.y);
            cout << fline;
            if(!DbleEqual(grid.z,0.))
            {	sprintf(fline," dz: %g",grid.z);
                cout << fline;
            }
            cout << endl;
            
            SetParticleLength(pointsPerCell);
            
            if(isAxisym)
                sprintf(fline,"Origin at R: %g Z: %g",xmin,ymin);
            else
                sprintf(fline,"Origin at x: %g y: %g",xmin,ymin);
            cout << fline;
            if(!DbleEqual(grid.z,0.))
            {	sprintf(fline," z: %g",zmin);
                cout << fline;
            }
            cout << endl;
        }
		
		if(DbleEqual(grid.z,0.) && !isAxisym)
		{	sprintf(fline,"Thickness: %g",zmin);
			cout << fline << endl;
		}

#ifdef _OPENMP
		if(isAxisym)
			sprintf(fline,"Patch Grid R: %d Z: %d",xpnum,ypnum);
		else
			sprintf(fline,"Patch Grid x: %d y: %d",xpnum,ypnum);
		cout << fline;
		if(!DbleEqual(grid.z,0.))
		{	sprintf(fline," z: %d",zpnum);
			cout << fline;
		}
		cout << endl;
#endif
	}
	else
		cout << "Non-orthogonal grid" << endl;
}

// output contact method by displacements or position with a cutoff
void MeshInfo::OutputContactByDisplacements(void)
{
    if(contactByDisplacements)
		cout << "   (normal cod from displacements)" << endl;
	else
		cout << "   (normal cod from position with contact when separated less than " << positionCutoff
                << " of cell)" << endl;
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

// For structured with equal element sizes only, find element from location and return result (1-based element number)
// Call code must be sure it is structured grid with equal element sizes
int MeshInfo::FindElementFromPoint(Vector *pt)
{
    int theElem;
    
    int col = (int)((pt->x-xmin)/grid.x);		// zero-based column # from 0 to horiz-1
	if(col<0 || col>=horiz)
	{	if(pt->x == xmin+horiz*grid.x)
			col = horiz-1;
		else
        {   char msg[100];
            sprintf(msg,"column for point (%lf,%lf,%lf)",pt->x,pt->y,pt->z);
			throw CommonException(msg,"");
        }
	}
    
    int row = (int)((pt->y-ymin)/grid.y);        // zero-based row # from 0 to vert-1
	if(row<0 || row>=vert)
	{	if(pt->y == ymin+vert*grid.y)
			row = vert-1;
		else
        {   char msg[100];
            sprintf(msg,"row for point (%lf,%lf,%lf)",pt->x,pt->y,pt->z);
            throw CommonException(msg,"");
        }
	}
    
    // 3D
    if(grid.z > 0.)
    {   int zrow = (int)((pt->z-zmin)/grid.z);   // zero-based row # from 0 to depth-1
		if(zrow<0 || zrow>=depth)
		{	if(pt->z == zmin+depth*grid.z)
				zrow = depth-1;
			else
            {   char msg[100];
                sprintf(msg,"rank for point (%lf,%lf,%lf)",pt->x,pt->y,pt->z);
                throw CommonException(msg,"");
            }
		}
        theElem = horiz*(zrow*vert + row) + col + 1;
    }
    
    // 2D
    else
    {   theElem = row*horiz + col + 1;
    }
    
    // return result
    return theElem;
}

// Create the patches for the grid
// Return pointer to a 0-based listed or patches[0] ... pathes[numProcs-1]
// Return NULL on memory error
GridPatch **MeshInfo::CreatePatches(int np,int numProcs)
{
	// serial or single thread is simpler
	if(numProcs<=1)
	{	return CreateOnePatch(np);
	}
    
    // custom patching
    if(fmobj->dflag[3]>0)
    {   // decode from flag as xxyyzz
        xpnum = fmobj->dflag[3]/10000;
        ypnum = (fmobj->dflag[3]-xpnum*10000)/100;
        zpnum = (fmobj->dflag[3]-xpnum*10000-ypnum*100);
        if(xpnum*ypnum*zpnum != numProcs) xpnum = -1;
    }
    else
        xpnum = -1;
	
	// get prime factors in ascending order
	vector<int> factors;
	unsigned ndim = np==THREED_MPM ? 3 : 2 ;
	PrimeFactors(numProcs,factors);
	while(factors.size()<ndim) factors.push_back(1);
	while(factors.size()>ndim)
	{	std::sort(factors.begin(),factors.end());
		int newFactor = factors[0]*factors[1];
		factors.erase(factors.begin());
		factors[0]=newFactor;
	}
	std::sort(factors.begin(),factors.end());
    
    if(xpnum<0)
    {   // convert to numbers of elements in each direction and use larger
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
    
    // get patch sizes
	xPatchSize = max(int(horiz/xpnum+.5),1);
	yPatchSize = max(int(vert/ypnum+.5),1);
    if(ndim==2)
        zPatchSize = 1;
    else
        zPatchSize = max(int(depth/zpnum+.5),1);
    
    // alloc space for patches - exit on memory error
    int totalPatches = xpnum*ypnum*zpnum;
    GridPatch **patch = (GridPatch **)malloc(totalPatches*sizeof(GridPatch));
    if(patch==NULL) return NULL;
	
	// create the patches
    int pnum=0;
	int i,j,k;
	int x1,x2,y1,y2,z1=1,z2;
	for(k=1;k<=zpnum;k++)
	{	z2 = k==zpnum ? depth : z1+zPatchSize-1;
		y1 = 1;
		for(j=1;j<=ypnum;j++)
		{	y2 = j==ypnum ? vert : y1+yPatchSize-1;
			x1 = 1;
			for(i=1;i<=xpnum;i++)
			{	x2 = i==xpnum ? horiz : x1+xPatchSize-1;
				
				// patch x1 to x2 and y1 to y2 (1 based)
				//cout << "\n- Patch " << pnum << ":" << i << "-" << j << "-" << k << ":";
                patch[pnum] = new GridPatch(x1,x2,y1,y2,z1,z2);
				if(!patch[pnum]->CreateGhostNodes()) return NULL;
                pnum++;
				
				x1 = x2+1;
			}
			y1 = y2+1;
		}
		z1 = z2+1;
	}
	
	// fill patches with particles
	int pn;
	for(int p=0;p<nmpms;p++)
	{	pn = GetPatchForElement(mpm[p]->ElemID());
		if(pn<0 || pn>=totalPatches)
		{	cout << "particle in element " << mpm[p]->ElemID() << " to patch " << pn << endl;
			continue;
		}
		patch[pn]->AddParticle(mpm[p]);
	}
	
    // return array of patches
    return patch;
}

// Create a single patch for the grid and patch has no ghost nodes
// Return pointer to a 0-based listed or patches[0]
// Return NULL on memory error
GridPatch **MeshInfo::CreateOnePatch(int np)
{
	// a single patch
	xpnum = ypnum = zpnum = 1;
	
	xPatchSize = horiz;
	yPatchSize = vert;
	zPatchSize = np==THREED_MPM ? depth : 1 ;
    
    // alloc space for patches - exit on memory error
    GridPatch **patch = (GridPatch **)malloc(sizeof(GridPatch));
    if(patch==NULL) return NULL;
	
	// one patch but no ghost nodes
	patch[0] = new GridPatch(1,xPatchSize,1,yPatchSize,1,np==THREED_MPM ? depth : 0);
	
	// fill patch with all particles
	for(int p=0;p<nmpms;p++) patch[0]->AddParticle(mpm[p]);
    
    // return array with the one patch
    return patch;
}

#pragma mark MeshInfo:Accessors

// get some info needed for symmetry BCs
double MeshInfo::GetParametersForBCs(int axis,double *gmin,double *gmax)
{
	if(axis==X_DIRECTION)
	{	*gmin = xmin;
		*gmax = xmax;
		return grid.x;
	}
	else if(axis==Y_DIRECTION)
	{	*gmin = ymin;
		*gmax = ymax;
		return grid.y;
	}
	else
	{	*gmin = zmin;
		*gmax = zmax;
		return grid.z;
	}
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
		int snum = iel % perSlice;		// number in slice (0 to horz*vert-1)
		col = snum % horiz;				// col 0 to horiz-1
		row = snum/horiz;				// zero based
		pcol = min(col/xPatchSize,xpnum-1);
		prow = min(row/yPatchSize,ypnum-1);
		prank = min(rank/zPatchSize,zpnum-1);
		return xpnum*(ypnum*prank + prow) + pcol;
	}
	
	// 2D
	col = iel % horiz;			// col 0 to horiz-1
	row = iel/horiz;			// zero based
	pcol = min(col/xPatchSize,xpnum-1);
	prow = min(row/yPatchSize,ypnum-1);
	return xpnum*prow + pcol;
}

// set grid style (zcell=0 if 2D grid)
// Options: style=false (NOT_CARTESIAN) means not aligned with x,y,z axes
//    style=VARIABLE_ORTHOGONAL_GRID
//      3D grid with unequal element sizes
//	  style=VARIABLE_RECTANGULAR_GRID
//		2D grid with unequal element sizes
//    style=SQUARE_GRID and all elements the same size as subsets, but adjust as follows
//      SQUARE_GRID = 2D and dx = dy
//      CUBIC_GRID = 3D and dx = dy = dz
//      RECTANGULAR_GRID = 2D and dx != dy
//      ORTHONGONAL_GRID = 3D and one of dx, dy, dz != others
//      set grid.i and diag.i for all these cases
void MeshInfo::SetCartesian(int style,double xcell,double ycell,double zcell)
{
	cartesian=style;
    if(style == VARIABLE_ORTHOGONAL_GRID || style == VARIABLE_RECTANGULAR_GRID)
    {   grid.z = zcell;
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
		
        // find normal vector on diagonal of the cell (currently not used)
		double diaglen=sqrt(grid.x*grid.x+grid.y*grid.y+grid.z*grid.z);
		diag.x=grid.x/diaglen;
		diag.y=grid.y/diaglen;
        diag.z=grid.z/diaglen;
	}
}

// set grid style (only set by <Grid> command). For 2D, d=0 and z is the specified thickness, or 1.0 by default,
//   or 1.0 if axisymmetric (so cell Volume is area)
// horiz is number of elements in x direction, vert is y direction, depth in z direction (or 0 in 2D)
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

// set grid length (as radius) per particle - requires cartesian grid, but need not be structured
void MeshInfo::SetParticleLength(int pointsPerCell)
{
	// surface length per particle on edge
	switch(pointsPerCell)
	{	case 1:
			lp = 1.0;
			break;
        case 9:
        case 27:
            // 9 is always 2D and 27 is always 3D
            lp = 1./3.;
            break;
        case 16:
            lp = 0.25;
            break;
        case 25:
            lp = 0.20;
            break;
		case 4:
		case 8:
		default:
			// 4 is always 2D and 8 is always 3D
			lp = 0.5;
			break;
	}
	
    // semi particle size in dimensioned uits
    part.x = lp*grid.x/2.;
	part.y = lp*grid.y/2.;
	part.z = lp*grid.z/2.;
}

// return cartesian setting
int MeshInfo::GetCartesian(void) { return cartesian; }

// see if 3D grid
bool MeshInfo::Is3DGrid(void) { return cartesian > BEGIN_3D_GRIDS; }

// return minimum cell dimension if the grid is cartesian (-1 if not)
// should find the minimum size in prelimnary calculations and store here
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

// horiz will be set >0 if generated by <Grid> command (but elements may have different sizes)
bool MeshInfo::IsStructuredGrid(void) { return horiz>0; }

// horiz will be set >0 if generated by <Grid> command and won't have variable element sizes
bool MeshInfo::IsStructuredEqualElementsGrid(void)
{	return horiz>0 && cartesian!=VARIABLE_ORTHOGONAL_GRID && cartesian!=VARIABLE_RECTANGULAR_GRID;
}

// number of points in the three directions (1 more than number of elements, ptz will be 1 if 2D)
// Feature that calls this method must require the problem to have a structured <Grid>
void MeshInfo::GetGridPoints(int *ptx,int *pty,int *ptz)
{	*ptx=horiz+1;
	*pty=vert+1;
	*ptz=depth+1;
}

// cell volume in mm^3
// Feature that calls this method must require the problem to have a structured <Grid>
double MeshInfo::GetCellVolume(void) { return cellVolume; }

// cell size in mm
// Feature that calls this method must require the problem to have a structured <Grid>
double MeshInfo::GetAverageCellSize(void) { return avgCellSize; }

// grid thickness (or -1. if not structured grid). For 3D returns z extent
// but not used in 3D
double MeshInfo::GetThickness(void)
{	if(horiz<=0) return -1.0;
	if(DbleEqual(grid.z,0.)) return zmin;
	return grid.z*depth;
}

// get thickness or 1.0 if not a structure grid
double MeshInfo::GetDefaultThickness()
{	double gthick=GetThickness();
	return gthick>0. ? gthick : 1.0 ;
}

// find cutoff distance for contact
double MeshInfo::GetNormalCODAdjust(Vector *norm,Vector *tang,double delt)
{
    // none if contact is by displacements
    if(contactByDisplacements) return 0.;
    
    // otherwise position cuttoff from fraction of perpendicular distance
    return positionCutoff * GetPerpendicularDistance(norm,tang,delt);
}

// find area of contact at a node
// delta is cod and norm is a normalized normal at the interface
// rawSurfaceArea is precalculated raw area
// tang, deln, and delt are outputs of these calculations (pointed must be provided)
double MeshInfo::InterfaceContactArea(Vector *delta,Vector *norm,double rawSurfaceArea,
									  Vector *tang,double *delnout,double *deltout)
{
	// perpendicular distance to correct contact area and contact by positions
    double dist;
    
    // normal displacement (norm is normalized) = delta . norm, subtract adjustment when using position
	// which have been precalculated
    double deln = DotVectors(delta,norm);
	
    // tangential vector in tang
    CopyVector(tang,delta);
    AddScaledVector(tang,norm,-deln);				// delta - deln (n) = dott (t)
    double delt=sqrt(DotVectors(tang,tang));
    if(!DbleEqual(delt,0.)) ScaleVector(tang,1/delt);
	
	// if using displacements, find dist  with tang
	if(GetContactByDisplacements())
	{	dist = mpmgrid.GetPerpendicularDistance(norm, tang, delt);
	}
	
	else
	{	// Cannot use tang here, which means 3D non regular will be less accurate
		// But, all regular and all 2D will be correct without tang
		// Future Goal: get hperp without needing tang for general 3D case
		// (for efficiency, call hperp method separately)
		dist = GetPerpendicularDistance(norm, NULL, 0.);
        deln -= positionCutoff*dist;
	}
 	
	*delnout = deln;
	*deltout = delt;
	return rawSurfaceArea/dist;
}

// find hperp distance used in contact calculations in interface force calculations
// Vector tang and magnitude delt only needed for 3D calculations. If not known
// or has zero magnitude, it will use normal vector method instead
double MeshInfo::GetPerpendicularDistance(Vector *norm,Vector *tang,double delt)
{
    double dist = grid.x;
    
    // Angled path correction method 1: hperp  is distance to ellipsoid through cell corners
    //    defined by tangent vector. In 3D, also multiply by distance to ellipsoid along
    //    n X t (which is along z axis for 2D)
    // In 2D and 3D the dist is equal to grid spacing if grid.x=grid.y=grid.z and therefore this
    //    whole block gets skipped
    // See JANOSU-6-60 and JANOSU-6-74 and method #1 in paper
    if(Is3DGrid())
    {   if(cartesian!=CUBIC_GRID)
        {   if(tang==NULL || DbleEqual(delt,0.))
            {   // rather then try to pick a tangent, use normal only
				// or use method #2 in paper and below
				double a=norm->x/mpmgrid.grid.x;
				double b=norm->y/mpmgrid.grid.y;
				double c=norm->z/mpmgrid.grid.z;
				dist = 1./sqrt(a*a + b*b + c*c);
            }
			else
			{	Vector t2;
				t2.x = norm->y*tang->z - norm->z*tang->y;
				t2.y = norm->z*tang->x - norm->x*tang->z;
				t2.z = norm->x*tang->y - norm->y*tang->x;
				double a1 = tang->x/grid.x;
				double b1 = tang->y/grid.y;
				double c1 = tang->z/grid.z;
				double a2 = t2.x/grid.x;
				double b2 = t2.y/grid.y;
				double c2 = t2.z/grid.z;
				dist = grid.x*grid.y*grid.z*sqrt((a1*a1 + b1*b1 + c1*c1)*(a2*a2 + b2*b2 + c2*c2));
			}
        }
    }
    else if(cartesian!=SQUARE_GRID)
    {   double a=grid.x*norm->x;
        double b=grid.y*norm->y;
        dist = sqrt(a*a + b*b);
    }

    // Angled path correction method 2: distance to ellipsoid along normal
    //      defined as hperp
    // See JANOSU-6-76 and method #2 is paper
	/*
    double a=norm->x/mpmgrid.grid.x;
    double b=norm->y/mpmgrid.grid.y;
    if(mpmgrid.Is3DGrid())
    {   double c=norm->z/mpmgrid.grid.z;
        dist = 1./sqrt(a*a + b*b + c*c);
    }
    else
        dist = 1./sqrt(a*a + b*b);
	*/

    // Angled path correction method 3 (in imperfect interface by cracks paper):
    //   Find perpendicular distance which gets smaller as interface tilts
    //   thus the effective surface area increases
    // See JANOSU-6-23 to 49 and method #3 in paper
	/*
    double a=fabs(mpmgrid.grid.x*norm->x);
    double b=fabs(mpmgrid.grid.y*norm->y);
    if(mpmgrid.Is3DGrid())
    {   // 3D has two cases
        double c=fabs(mpmgrid.grid.z*norm->z);
        dist = fmax(a,fmax(b,c));
        if(2.*dist < a+b+c)
        {   // need alternate formula in this case (i.e., Max(a,b,c) < sum of other two)
            dist = (1./4.)*(2./a + 2./b + 2/c - a/(b*c) - b/(a*c) - c/(a*b));
            dist = 1./dist;
        }
    }
    else
    {   // 2D just take maximum
        dist = fmax(a,b);
    }
	*/
	
    return dist;
}
	
// return current setting for contact method
bool MeshInfo::GetContactByDisplacements(void) { return contactByDisplacements; }
void MeshInfo::SetContactByDisplacements(bool newContact) { contactByDisplacements=newContact; }

// Cell size for grid with constant element size only (not checked here)
Vector MeshInfo::GetCellSize(void) { return grid; }
double MeshInfo::GetCellXSize(void) { return grid.x; }

// Particle size for grid with constant element size only
Vector MeshInfo::GetParticleSize(void) { return part; }
double MeshInfo::GetParticleXSize(void) { return part.x; }
double MeshInfo::GetParticleYSize(void) { return part.y; }
double MeshInfo::GetParticleZSize(void) { return part.z; }

// particle semi length in natural coordinates
double MeshInfo::GetParticleSemiLength(void) { return lp; }




