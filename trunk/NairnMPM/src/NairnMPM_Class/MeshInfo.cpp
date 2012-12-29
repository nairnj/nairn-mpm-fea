/*********************************************************************
    MeshInfo.cpp
    Nairn Research Group MPM Code
    
    Created by John Nairn on 5/16/06.
    Copyright (c) 2006, All rights reserved.
	
	Dependencies
		none
*********************************************************************/

#include "MeshInfo.hpp"

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
void MeshInfo::Output(int pointsPerCell)
{
    char fline[200];

	if(cartesian>0)
	{	if(horiz>0)
			cout << "Structured";
		else
			cout << "Unstructured";
    	sprintf(fline," orthogonal grid with dx: %g dy: %g",gridx,gridy);
		cout << fline;
		if(!DbleEqual(gridz,0.))
		{	sprintf(fline," dz: %g",gridz);
			cout << fline;
		}
		cout << endl;
		SetParticleLength(pointsPerCell);
    	sprintf(fline,"Origin at x: %g y: %g",xmin,ymin);
		cout << fline;
		if(!DbleEqual(gridz,0.))
		{	sprintf(fline," z: %g",zmin);
			cout << fline;
		}
		cout << endl;
		if(DbleEqual(gridz,0.))
		{	sprintf(fline,"Thickness: %g",zmin);
			cout << fline << endl;
		}
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
	int i=0;
	int col=num % horiz;
	int below=num-horiz;
	int above=num+horiz;
	
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
// NOT WRITTEN YET
void MeshInfo::ListOfNeighbors3D(int num,int *neighbor)
{
	if(horiz<=0)
	{	neighbor[0]=0;
		return;
	}
	
	int i=0;
	
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

#pragma mark MeshInfo:Accessors

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
	}
	else
	{	totalElems=horiz*vert;
		cellVolume=gridx*gridy*zmin;
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
// Vector tang and magnitude delt only needed for 3D calculations
double MeshInfo::GetPerpendicularDistance(Vector *norm,Vector *tang,double delt)
{
    double dist = gridx;
    
    // Angled path correction method 1: hperp  is distance to ellipsoid through cell corners
    //    defined by tangent vector. In 3D, also multiply by distance to ellipsoid along
    //    n X t (which is along z axis for 2D)
    // In 2D and 3D the dist is equal to grid spacing if gridx=gridy=gridz and therefore this
    //    whole block gets skipped
    // See JANOSU-6-60 and JANOSU-6-74
    if(Is3DGrid())
    {   if(cartesian!=CUBIC_GRID)
        {   if(DbleEqual(delt,0.))
            {   // pick any tangent vector
                tang->z = 0.;
                if(!DbleEqual(norm->x,0.0) || !DbleEqual(norm->y,0.0))
                {   tang->x = norm->y;
                    tang->y = -norm->x;
                }
                else
                {   // norm = (0,0,1)
                    tang->x = 1.;
                    tang->y = 0.;
                }
            }
            Vector t2;
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
    else if(cartesian!=SQUARE_GRID)
    {   double a=gridx*norm->x;
        double b=gridy*norm->y;
        dist = sqrt(a*a + b*b);
    }

    // Angled path correction method 2: distance to ellipsoid along normal
    //      defined as hperp
    // See JANOSU-6-76
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
    // See JANOSU-6-23 to 49
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


	

