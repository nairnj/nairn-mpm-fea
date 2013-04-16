/*********************************************************************
    MeshInfo.cpp
    Nairn Research Group MPM Code
    
    Created by John Nairn on 5/16/06.
    Copyright (c) 2006, All rights reserved.
	
	Dependencies
		none
*********************************************************************/

#include "NairnMPM_Class/MeshInfo.hpp"

// NEWINCLUDE
#include "Patches/GridPatch.hpp"
#include "MPM_Classes/MPMBase.hpp"

// global class for grid information
MeshInfo mpmgrid;

#pragma mark MeshInfo:Constructors and Destructors

MeshInfo::MeshInfo(void)
{
	cartesian=UNKNOWN_GRID;
	horiz=0;                        // also flag that used <Grid> command (i.e. structured grid)
	contactByDisplacements=TRUE;	// contact by displacements
	positionCutoff=0.8;             // element fraction when contact by positions
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
		
		if(isAxisym)
			sprintf(fline," orthogonal grid with dR: %g dZ: %g",gridx,gridy);
		else
			sprintf(fline," orthogonal grid with dx: %g dy: %g",gridx,gridy);
		cout << fline;
		if(!DbleEqual(gridz,0.))
		{	sprintf(fline," dz: %g",gridz);
			cout << fline;
		}
		cout << endl;
		
		SetParticleLength(pointsPerCell);
		
		if(isAxisym)
			sprintf(fline,"Origin at R: %g Z: %g",xmin,ymin);
		else
			sprintf(fline,"Origin at x: %g y: %g",xmin,ymin);
		cout << fline;
		if(!DbleEqual(gridz,0.))
		{	sprintf(fline," z: %g",zmin);
			cout << fline;
		}
		cout << endl;
		
		if(DbleEqual(gridz,0.) && !isAxisym)
		{	sprintf(fline,"Thickness: %g",zmin);
			cout << fline << endl;
		}

#ifdef _OPENMP
		if(isAxisym)
			sprintf(fline,"Patch Grid R: %d Z: %d",xpnum,ypnum);
		else
			sprintf(fline,"Patch Grid x: %d y: %d",xpnum,ypnum);
		cout << fline;
		if(!DbleEqual(gridz,0.))
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

// If have a structured grid, get the 27 (3D) neighbor elements. The 1-based
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
    // 1-based element number below, at, and above element num
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
			if(j!=num) neighbor[i++]=j;
			neighbor[i++]=j-1;
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
			if(j!=num) neighbor[i++]=j;
			neighbor[i++]=j-1;
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

// For structured grid only, find element from location and return result (1-based element number)
// throws empty exception if not a structured grid (horiz<0) or if point is not in the grid
int MeshInfo::FindElementFromPoint(Vector *pt)
{
    int theElem = 0;
    
    // error if not structured grid
    if(horiz<0) throw "";
    
    int col = (int)((pt->x-xmin)/gridx);		// zero-based column # from 0 to horiz-1
	if(col<0 || col>=horiz)
	{	if(pt->x == xmin+horiz*gridx)
			col = horiz-1;
		else
			throw "";
	}
    int row = (int)((pt->y-ymin)/gridy);
	if(row<0 || row>=vert)
	{	if(pt->y == ymin+vert*gridy)
			row = vert-1;
		else
			throw "";
	}
    
    // 3D
    if(gridz > 0.)
    {   int zrow = (int)((pt->z-zmin)/gridz);
		if(zrow<0 || zrow>=depth)
		{	if(pt->z == zmin+depth*gridz)
				zrow = depth-1;
			else
				throw "";
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
	
	// convert to numbers of elements in each direction and use larger
	// prime factors in the longer directions
	if(ndim==2)
	{	if(horiz>vert)
		{	xpnum = factors[1];
			ypnum = factors[0];
		}
		else
		{	xpnum = factors[0];
			ypnum = factors[1];
		}
		zpnum=1;
		zPatchSize=1;
	}
	else
	{	if(horiz>vert && horiz>depth)
		{	xpnum = factors[2];
			if(vert>depth)
			{	ypnum = factors[1];
				zpnum = factors[0];
			}
			else
			{	zpnum = factors[0];
				ypnum = factors[1];
			}
		}
		else if(vert>horiz && vert>depth)
		{	ypnum = factors[2];
			if(horiz>depth)
			{	xpnum = factors[1];
				zpnum = factors[0];
			}
			else
			{	zpnum = factors[0];
				xpnum = factors[1];
			}
		}
		else
		{	zpnum = factors[2];
			if(horiz>vert)
			{	xpnum = factors[1];
				ypnum = factors[0];
			}
			else
			{	ypnum = factors[0];
				xpnum = factors[1];
			}
		}
		zPatchSize = max(int(depth/zpnum+.5),1);
	}
	xPatchSize = max(int(horiz/xpnum+.5),1);
	yPatchSize = max(int(vert/ypnum+.5),1);
    
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
				cout << "\n- Patch " << pnum << ":" << i << "-" << j << "-" << k << ":";
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

// Create a single patches for the grid and patch has no ghost nodes
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
	
	// one patch an no ghost nodes
	patch[0] = new GridPatch(1,xPatchSize,1,yPatchSize,1,zPatchSize);
	if(!patch[0]->CreateGhostNodes()) return NULL;
	
	// fill patch with particles
	for(int p=0;p<nmpms;p++) patch[0]->AddParticle(mpm[p]);
    
    // return array with the one patch
    return patch;
}

#pragma mark MeshInfo:Accessors

// given zero based element number, return patch number
int MeshInfo::GetPatchForElement(int iel)
{
	int row,col,rank;
	int prow,pcol,prank;
	
	if(depth>0)
	{	// 3D
		int perSlice = horiz*vert;		// number in each slice
		rank = iel/perSlice;
		int snum = iel % perSlice;		// number in slice (0 to horz*vert-1)
		col = snum % horiz;				// col 0 to horiz-1
		row = snum/horiz;				// zero based
		pcol = min(col/xPatchSize,xpnum-1);
		prow = min(row/yPatchSize,ypnum-1);
		prank = min(rank/zPatchSize,zpnum-1);
		return xpnum*ypnum*prank + xpnum*prow + pcol;
	}
	
	// 2D
	col = iel % horiz;			// col 0 to horiz-1
	row = iel/horiz;			// zero based
	pcol = min(col/xPatchSize,xpnum-1);
	prow = min(row/yPatchSize,ypnum-1);
	return xpnum*prow + pcol;
}

// set grid style (zcell=0 if 2D grid)
// Options: style=0 (NOT_CARTESIAN) means not aligned with x,y,z axes
//    SQUARE_GRID = 2D and dx = dy
//    CUBIC_GRID = 3D and dx = dy = dz
//    RECTANGULAR_GRID = 2D and dx != dy
//    ORTHONGONAL_GRID = 3D and one of dx, dy, dz != others
// If cartensian, set gridi and celli
void MeshInfo::SetCartesian(int style,double xcell,double ycell,double zcell)
{
	cartesian=style;
	if(style>0)
	{	gridx=xcell;
		gridy=ycell;
		gridz=zcell;
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
		
        // find normal vector on diagonal of the cell
		double diag=sqrt(gridx*gridx+gridy*gridy+gridz*gridz);
		diagx=gridx/diag;
		diagy=gridy/diag;
        diagz=gridz/diag;
	}
}

// set grid style (only set by <Grid> command). For 2D, d=0 and z is the specified thickness, or 1.0 by default,
//   or 1.0 if axisymmetric (so cell Volume is area)
// horiz is number of elements in x direction, vert is y direction, depth in z direction (or 0 in 2D)
void MeshInfo::SetElements(int h,int v,int d,double x,double y,double z)
{	horiz=h;
	vert=v;
	depth=d;
	xmin=x;
	ymin=y;
	zmin=z;
	
	// spacing between nodes numbers along each axis
	xplane=1;
	yplane=horiz+1;
	zplane=(horiz+1)*(vert+1);
	if(depth>0)
	{	totalElems=horiz*vert*depth;
		cellVolume=gridx*gridy*gridz;
        avgCellSize=(gridx+gridy+gridz)/3.;
	}
	else
	{	totalElems=horiz*vert;
		cellVolume=gridx*gridy*zmin;
        avgCellSize=(gridx+gridy)/2.;
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
	
	partx = lp*gridx/2.;
	party = lp*gridy/2.;
	partz = lp*gridz/2.;
}

// return cartesian setting
int MeshInfo::GetCartesian(void) { return cartesian; }

// see if 3D grid
bool MeshInfo::Is3DGrid(void) { return cartesian > BEGIN_3D_GRIDS; }

// return minimum cell dimension if the grid is cartesian (-1 if not)
double MeshInfo::GetMinCellDimension(void)
{
	double minSize;
	
	switch(cartesian)
	{	case SQUARE_GRID:
		case CUBIC_GRID:
			minSize=gridx;
			break;
			
		case RECTANGULAR_GRID:
			minSize=fmin(gridx,gridy);
			break;
			
		case ORTHOGONAL_GRID:
			minSize=fmin(gridx,gridy);
			minSize=fmin(minSize,gridz);
			break;
			
		default:
			minSize=-1.;
			break;
	}
	
	return minSize;
}

// can mesh do GIMP - must be Cartesian and structured using <Grid> command
int MeshInfo::CanDoGIMP(void)
{
	// must have used <Grid> to set up mesh
	if(cartesian<=0) return FALSE;
	if(horiz==0) return FALSE;
	return TRUE;
}

// horiz will be set >0 if generated by <Grid> command
bool MeshInfo::IsStructuredGrid(void) { return horiz>0; }

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

// grid thickness (or -1. if not structure grid). For 3D return z extent
double MeshInfo::GetThickness(void)
{	if(horiz<=0) return -1.0;
	if(DbleEqual(gridz,0.)) return zmin;
	return gridz*depth;
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

// find hperp distance used in contact calculations in interface force calculations
// Vector tang and magnitude delt only needed for 3D calculations. If not known
// or has zero magnitude, it will use normal vector method instead
double MeshInfo::GetPerpendicularDistance(Vector *norm,Vector *tang,double delt)
{
    double dist = gridx;
    
    // Angled path correction method 1: hperp  is distance to ellipsoid through cell corners
    //    defined by tangent vector. In 3D, also multiply by distance to ellipsoid along
    //    n X t (which is along z axis for 2D)
    // In 2D and 3D the dist is equal to grid spacing if gridx=gridy=gridz and therefore this
    //    whole block gets skipped
    // See JANOSU-6-60 and JANOSU-6-74 and method #1 in paper
    if(Is3DGrid())
    {   if(cartesian!=CUBIC_GRID)
        {   if(tang==NULL || DbleEqual(delt,0.))
            {   // rather then try to pick a tangent, use normal only
				// or use method #2 in paper and below
				double a=norm->x/mpmgrid.gridx;
				double b=norm->y/mpmgrid.gridy;
				double c=norm->z/mpmgrid.gridz;
				dist = 1./sqrt(a*a + b*b + c*c);
            }
			else
			{	Vector t2;
				t2.x = norm->y*tang->z - norm->z*tang->y;
				t2.y = norm->z*tang->x - norm->x*tang->z;
				t2.z = norm->x*tang->y - norm->y*tang->x;
				double a1 = tang->x/gridx;
				double b1 = tang->y/gridy;
				double c1 = tang->z/gridz;
				double a2 = t2.x/gridx;
				double b2 = t2.y/gridy;
				double c2 = t2.z/gridz;
				dist = gridx*gridy*gridz*sqrt((a1*a1 + b1*b1 + c1*c1)*(a2*a2 + b2*b2 + c2*c2));
			}
        }
    }
    else if(cartesian!=SQUARE_GRID)
    {   double a=gridx*norm->x;
        double b=gridy*norm->y;
        dist = sqrt(a*a + b*b);
    }

    // Angled path correction method 2: distance to ellipsoid along normal
    //      defined as hperp
    // See JANOSU-6-76 and method #2 is paper
	/*
    double a=norm->x/mpmgrid.gridx;
    double b=norm->y/mpmgrid.gridy;
    if(mpmgrid.Is3DGrid())
    {   double c=norm->z/mpmgrid.gridz;
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
    double a=fabs(mpmgrid.gridx*norm->x);
    double b=fabs(mpmgrid.gridy*norm->y);
    if(mpmgrid.Is3DGrid())
    {   // 3D has two cases
        double c=fabs(mpmgrid.gridz*norm->z);
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


	

