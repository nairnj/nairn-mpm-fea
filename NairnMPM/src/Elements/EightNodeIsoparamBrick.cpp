/******************************************************************
	EightNodeIsoparamBrick.cpp
	nairn-mpm-fea

	Created by John Nairn on 7/20/06.
	Copyright 2006 RSAC Software. All rights reserved.
 ******************************************************************/

#include "stdafx.h"
#include "Elements/EightNodeIsoparamBrick.hpp"
#include "Nodes/NodalPoint.hpp"
#ifdef MPM_CODE
	#include "NairnMPM_Class/MeshInfo.hpp"
#endif

// Local globals
static double xii[8]={-1.,1.,1.,-1.,-1.,1.,1.,-1.};
static double eti[8]={-1.,-1.,1.,1.,-1.,-1.,1.,1.};
static double zti[8]={-1.,-1.,-1.,-1.,1.,1.,1.,1.};

// globals for GIMP node locations
#ifdef MPM_CODE
static double g3xii[64]={-1.,1.,1.,-1.,-1.,1.,1.,-1.,
						-3.,-1.,1.,3.,3.,3.,3.,1.,-1.,-3.,-3.,-3.,
						-3.,-1.,1.,3.,3.,3.,3.,1.,-1.,-3.,-3.,-3.,
						-1.,1.,1.,-1.,-3.,-1.,1.,3.,3.,3.,3.,1.,-1.,-3.,-3.,-3.,
						-1.,1.,1.,-1.,-3.,-1.,1.,3.,3.,3.,3.,1.,-1.,-3.,-3.,-3.};
static double g3eti[64]={-1.,-1.,1.,1.,-1.,-1.,1.,1.,
						-3.,-3.,-3.,-3.,-1.,1.,3.,3.,3.,3.,1.,-1.,
						-3.,-3.,-3.,-3.,-1.,1.,3.,3.,3.,3.,1.,-1.,
						-1.,-1.,1.,1.,-3.,-3.,-3.,-3.,-1.,1.,3.,3.,3.,3.,1.,-1.,
						-1.,-1.,1.,1.,-3.,-3.,-3.,-3.,-1.,1.,3.,3.,3.,3.,1.,-1.};
static double g3zti[64]={-1.,-1.,-1.,-1.,1.,1.,1.,1.,
						-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,
						1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
						-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,
						3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.};
// = int(0.5*(g3xii[i]+1)), int(0.5*(g3etai[i]+1)), and int(0.5*(g3zti[i]+1))
static int x3off[64]={0,1,1,0,0,1,1,0,
					-1,0,1,2,2,2,2,1,0,-1,-1,-1,
					-1,0,1,2,2,2,2,1,0,-1,-1,-1,
					0,1,1,0,-1,0,1,2,2,2,2,1,0,-1,-1,-1,
					0,1,1,0,-1,0,1,2,2,2,2,1,0,-1,-1,-1};
static int y3off[64]={0,0,1,1,0,0,1,1,
					-1,-1,-1,-1,0,1,2,2,2,2,1,0,
					-1,-1,-1,-1,0,1,2,2,2,2,1,0,
					0,0,1,1,-1,-1,-1,-1,0,1,2,2,2,2,1,0,
					0,0,1,1,-1,-1,-1,-1,0,1,2,2,2,2,1,0};
static int z3off[64]={0,0,0,0,1,1,1,1,
					0,0,0,0,0,0,0,0,0,0,0,0,
					1,1,1,1,1,1,1,1,1,1,1,1,
					-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
					2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
#endif

#pragma mark EightNodeIsoparamBrick::Constructors and Destructor

// main constructor - pass up the chain
EightNodeIsoparamBrick::EightNodeIsoparamBrick(int eNum,int *eNode) : ElementBase3D(eNum,eNode)
{
	nodes[8]=eNode[0];
}

#pragma mark EightNodeIsoparamBrick::Methods

// Get shape functions, but does note get derivatives. All calls for derivative
// must use second method which assume orthogonal elements
void EightNodeIsoparamBrick::ShapeFunction(Vector *xi,int getDeriv,
		double *sfxn,double *xiDeriv,double *etaDeriv,Vector *eNodes,
                double *outDetjac,double *outAsr,double *asbe) const
{
    double temp1,temp2,temp3;
    int i;
    
    // shape function
    for(i=0;i<8;i++)
    {	temp1=(1.+xii[i]*xi->x);
        temp2=(1.+eti[i]*xi->y);
        temp3=(1.+zti[i]*xi->z);
        sfxn[i]=temp1*temp2*temp3/8.;
    }
}

// get just shape functions and optionally derivative wrt x and y
// General for shape funciton, but derivatives assumes a rectangular mesh
void EightNodeIsoparamBrick::ShapeFunction(Vector *xi,int getDeriv,double *sfxn,
											double *xDeriv,double *yDeriv,double *zDeriv) const
{
    double temp1,temp2,temp3;
    int i;
    
    // shape function
    for(i=0;i<8;i++)
    {	temp1=(1.+xii[i]*xi->x);
        temp2=(1.+eti[i]*xi->y);
        temp3=(1.+zti[i]*xi->z);
        sfxn[i]=0.125*temp1*temp2*temp3;
		if(getDeriv)
		{	xDeriv[i]=0.25*xii[i]*temp2*temp3/GetDeltaX();
			yDeriv[i]=0.25*eti[i]*temp1*temp3/GetDeltaY();
			zDeriv[i]=0.25*zti[i]*temp1*temp2/GetDeltaZ();
		}
	}
}

// get B2SPLINE shape functions and optionally derivatives wrt x and y, but derivatives only work
// if it is a regular array. Shape functions are general
// For axisymmetric MPM, make sure zDeriv is not NULL and load with shape function/rp
void EightNodeIsoparamBrick::SplineShapeFunction(int *nds,Vector *xi,int getDeriv,double *sfxn,
										   double *xDeriv,double *yDeriv,double *zDeriv) const
{
	// constants
	double inv_dx = 0.,inv_dy=0.,inv_dz=0.;
	if(getDeriv)
	{   inv_dx = 1./GetDeltaX();
		inv_dy = 1./GetDeltaY();
		inv_dz = 1./GetDeltaZ();
	}
	
	// get cell shape functions
	double sx,sy,sz,arg,etai,netai,zetai,temp1,temp2,temp3;
	int i=0;
	for(int id=0;id<64;id++)
	{	// xi direction
		etai = xi->x - g3xii[id];
		temp1=fabs(etai);
		if(temp1>=3.) continue;
		
		// y direction
		netai = xi->y - g3eti[id];
		temp2=fabs(netai);
		if(temp2>=3.) continue;
		
		// z direction
		zetai = xi->z - g3zti[id];
		temp3=fabs(zetai);
		if(temp3>=3.) continue;
		
		// find shape function i
		
		// x direction
		if(temp1<=1.0)
			sx = 0.25*(3.-etai*etai);
		else
		{	arg = 3.-temp1;
			sx = 0.125*arg*arg;
		}
		
		// y direction
		if(temp2<=1.0)
			sy = 0.25*(3.-netai*netai);
		else
		{	arg = 3.-temp2;
			sy = 0.125*arg*arg;
		}
		
		// z direction
		if(temp3<=1.0)
			sz = 0.25*(3.-zetai*zetai);
		else
		{	arg = 3.-temp3;
			sz = 0.125*arg*arg;
		}

		// combine
		sfxn[i] = sx*sy*sz;
		
		// Note that derivative is (2/dx)dN/detai and (2/dy)dN/dnetai
		if(getDeriv)
		{	// x direction
			if(temp1<=1.0)
				xDeriv[i] = -etai*sy*sz*inv_dx;
			else if(etai>=0.)
				xDeriv[i] = 0.5*(etai-3)*sy*sz*inv_dx;
			else
				xDeriv[i] = 0.5*(etai+3)*sy*sz*inv_dx;
			
			// y direction
			if(temp2<=1.0)
				yDeriv[i] = -netai*sx*sz*inv_dy;
			else if(netai>=0.)
				yDeriv[i] = 0.5*(netai-3)*sx*sz*inv_dy;
			else
				yDeriv[i] = 0.5*(netai+3)*sx*sz*inv_dy;
			
			// z direction
			if(temp3<=1.0)
				zDeriv[i] = -zetai*sx*sy*inv_dz;
			else if(zetai>=0.)
				zDeriv[i] = 0.5*(zetai-3)*sx*sy*inv_dz;
			else
				zDeriv[i] = 0.5*(zetai+3)*sx*sy*inv_dz;
		}
		
		// assign node number
		i++;
		nds[i] = nodes[0] + x3off[id]*mpmgrid.xplane + y3off[id]*mpmgrid.yplane + z3off[id]*mpmgrid.zplane;
	}
	
	// set number
	// this better be <=27 or bad things happen
	nds[0] = i;
}

// see if point is in this element (assumes rectangular)
short EightNodeIsoparamBrick::PtInElement(Vector &pt) const
{	if(pt.x<xmin || pt.x>=xmax) return FALSE;
	if(pt.y<ymin || pt.y>=ymax) return FALSE;
	if(pt.z<zmin || pt.z>=zmax) return FALSE;
	return TRUE;
}

// see if this element is cube in cartesion coordinates and return TRUE or FALSE
// if cube, set dx, dy, dz to element dimensions
int EightNodeIsoparamBrick::Orthogonal(double *dx,double *dy,double *dz)
{
	int i;
	double xdel,ydel,zdel;
	
	*dx=0.;
	*dy=0.;
	*dz=0.;
	
	// front face
	for(i=0;i<3;i++)
    {	xdel=fabs(nd[nodes[i+1]]->x-nd[nodes[i]]->x);
		ydel=fabs(nd[nodes[i+1]]->y-nd[nodes[i]]->y);
		zdel=fabs(nd[nodes[i+1]]->z-nd[nodes[i]]->z);
		if(xdel>=1.e-12 && ydel>=1.e-12 && zdel>=1.e-12) return FALSE;
		*dx=fmax(xdel,*dx);
		*dy=fmax(ydel,*dy);
		*dz=fmax(zdel,*dz);
	}
    xdel=fabs(nd[nodes[0]]->x-nd[nodes[3]]->x);
	ydel=fabs(nd[nodes[0]]->y-nd[nodes[3]]->y);
	zdel=fabs(nd[nodes[0]]->z-nd[nodes[3]]->z);
	if(xdel>=1.e-12 && ydel>=1.e-12 && zdel>=1.e-12) return FALSE;
	*dx=fmax(xdel,*dx);
	*dy=fmax(ydel,*dy);
	*dz=fmax(zdel,*dz);
	
	// back face
	for(i=4;i<7;i++)
    {	xdel=fabs(nd[nodes[i+1]]->x-nd[nodes[i]]->x);
		ydel=fabs(nd[nodes[i+1]]->y-nd[nodes[i]]->y);
		zdel=fabs(nd[nodes[i+1]]->z-nd[nodes[i]]->z);
		if(xdel>=1.e-12 && ydel>=1.e-12 && zdel>=1.e-12) return FALSE;
		*dx=fmax(xdel,*dx);
		*dy=fmax(ydel,*dy);
		*dz=fmax(zdel,*dz);
	}
    xdel=fabs(nd[nodes[4]]->x-nd[nodes[7]]->x);
	ydel=fabs(nd[nodes[4]]->y-nd[nodes[7]]->y);
	zdel=fabs(nd[nodes[4]]->z-nd[nodes[7]]->z);
	if(xdel>=1.e-12 && ydel>=1.e-12 && zdel>=1.e-12) return FALSE;
	*dx=fmax(xdel,*dx);
	*dy=fmax(ydel,*dy);
	*dz=fmax(zdel,*dz);
	
	// side edges
	for(i=0;i<=3;i++)
    {	xdel=fabs(nd[nodes[i+4]]->x-nd[nodes[i]]->x);
		ydel=fabs(nd[nodes[i+1]]->y-nd[nodes[i]]->y);
		zdel=fabs(nd[nodes[i+1]]->z-nd[nodes[i]]->z);
		if(xdel>=1.e-12 && ydel>=1.e-12 && zdel>=1.e-12) return FALSE;
		*dx=fmax(xdel,*dx);
		*dy=fmax(ydel,*dy);
		*dz=fmax(zdel,*dz);
	}
	
	return TRUE;
}

// find dimensionless position, but assumes an orthogonal element
void EightNodeIsoparamBrick::GetXiPos(const Vector *pos,Vector *xipos) const
{
	xipos->x=(2.*pos->x-xmin-xmax)/GetDeltaX();
	xipos->y=(2.*pos->y-ymin-ymax)/GetDeltaY();
	xipos->z=(2.*pos->z-zmin-zmax)/GetDeltaZ();
}

#pragma mark EightNodeIsoparamBrick::GIMP Methods

// get GIMP shape functions and optionally derivatives wrt x and y
// assumed to be properly numbered regular 3D array
// input *xi position in element coordinate
// output number of nodes in nds[0] and node numbers in nds[1] to nds[nds[0]]
void EightNodeIsoparamBrick::GimpShapeFunction(Vector *xi,int *nds,int getDeriv,double *sfxn,
											   double *xDeriv,double *yDeriv,double *zDeriv,Vector &lp) const
{
	double xp,yp,zp,Svpx,Svpy,Svpz,dSvpx,dSvpy,dSvpz,xsign,ysign,zsign,argx=0.,argy=0.,argz=0.;
	
	// L is the cell spacing, 2*lpi is the current particle size (dimensionless for range -1 to 1).
	// The deformation of the particle is not considered yet.
    double q1x = 2.-lp.x, q2x = 2.+lp.x;
    double q1y = 2.-lp.y, q2y = 2.+lp.y;
    double q1z = 2.-lp.z, q2z = 2.+lp.z;

	// Pre-compute expensive divisions
	double inv_size_x = 1. / (4.*lp.x);
	double inv_size_y = 1. / (4.*lp.y);
	double inv_size_z = 1. / (4.*lp.y);
	double inv_dx = 0;
	double inv_dy = 0;
	double inv_dz = 0;
	if (getDeriv) {
		inv_dx = 2.0 / GetDeltaX();
		inv_dy = 2.0 / GetDeltaY();
		inv_dz = 2.0 / GetDeltaZ();
	}
	
	int i=0;
	for(int id=0;id<64;id++)
	{	// x direction
		xp=fabs(xi->x-g3xii[id]);
		if(xp>=q2x) continue;
		
		// y direction
		yp=fabs(xi->y-g3eti[id]);
		if(yp>=q2y) continue;
		
		// z direction
		zp=fabs(xi->z-g3zti[id]);
		if(zp>=q2z) continue;
		
		// shape function i
		
		if(xp<lp.x)
			Svpx = ((4.-lp.x)*lp.x-xp*xp)*inv_size_x;	// if lp=0.5: -(4.*xp*xp-7.)/8.;
		else if(xp<=q1x)
			Svpx = 0.5*(2.-xp);
		else
		{	argx = (q2x-xp)*inv_size_x;			// if lp=0.5: (5.-2.*xp)/4
			Svpx = 2.*lp.x*argx*argx;				// if lp=0.5: (5.-2.*xp)^2/16
		}
		
		if(yp<lp.y)
			Svpy = ((4.-lp.y)*lp.y-yp*yp)*inv_size_y;	// if lp=0.5: -(4.*yp*yp-7.)/8.;
		else if(yp<=q1y)
			Svpy = 0.5*(2.-yp);
		else
		{	argy = (q2y-yp)*inv_size_y;			// if lp=0.5: (5.-2.*yp)/4
			Svpy = 2.*lp.y*argy*argy;				// if lp=0.5: (5.-2.*yp)^2/16
		}
		
		if(zp<lp.z)
			Svpz = ((4.-lp.z)*lp.z-zp*zp)*inv_size_z;	// if lp=0.5: -(4.*zp*zp-7.)/8.;
		else if(zp<=q1z)
			Svpz = 0.5*(2.-zp);
		else
		{	argz = (q2z-zp)*inv_size_z;			// if lp=0.5: (5.-2.*zp)/4
			Svpz = 2.*lp.z*argz*argz;				// if lp=0.5: (5.-2.*zp)^2/16
		}
		
		sfxn[i] = Svpx*Svpy*Svpz;
		
		// find shape function gradient at (xp,yp,zp)
		if(getDeriv)
		{	xsign = xi->x>g3xii[id] ? 1. : -1.;
			ysign = xi->y>g3eti[id] ? 1. : -1.;
			zsign = xi->z>g3zti[id] ? 1. : -1.;
			
			if(xp<lp.x)
				dSvpx = -xp/(2.*lp.x);			// if lp=0.5: -xp
			else if(xp<=q1x)
				dSvpx = -0.5;
			else
				dSvpx = -argx;
			
			if(yp<lp.y)
				dSvpy = -yp/(2.*lp.y);			// if lp=0.5: -yp
			else if(yp<=q1y)
				dSvpy = -0.5;
			else
				dSvpy = -argy;
			
			if(zp<lp.z)
				dSvpz = -zp/(2.*lp.z);			// if lp=0.5: -zp;
			else if(zp<=q1z)
				dSvpz = -0.5;
			else
				dSvpz = -argz;
			
			xDeriv[i] = xsign*dSvpx*Svpy*Svpz*inv_dx;
			yDeriv[i] = ysign*Svpx*dSvpy*Svpz*inv_dy;
			zDeriv[i] = zsign*Svpx*Svpy*dSvpz*inv_dz;
		}
		
		// set node (and better have room)
		i++;
		nds[i] = nodes[0] + x3off[id]*mpmgrid.xplane + y3off[id]*mpmgrid.yplane + z3off[id]*mpmgrid.zplane;
	}
	
	// number of nodes (better be less than space available in nds[])
	nds[0] = i;
}

// get GIMP shape functions and optionally derivatives wrt x and y
// assumed to be properly numbered regular 3D array
// input *xi position in element coordinate and ndIDs[0]... is which nodes (0-63)
void EightNodeIsoparamBrick::BGimpShapeFunction(Vector *xi,int *nds,int getDeriv,double *sfxn,
											   double *xDeriv,double *yDeriv,double *zDeriv,Vector &lp) const
{
	double xp,yp,zp,Svpx,Svpy,Svpz,dSvpx,dSvpy,dSvpz,xsign,ysign,zsign,arg;
	
	// L is the cell spacing, 2*lpi is the current particle size (dimensionless range -1 to 1).
	// Breakpoints on positive side of the node
	double b1x = 1.-lp.x,b2x = 1.+lp.x,b3x = 3.-lp.x,b4x = 3.+lp.x;
	double b1y = 1.-lp.y,b2y = 1.+lp.y,b3y = 3.-lp.y,b4y = 3.+lp.y;
	double b1z = 1.-lp.z,b2z = 1.+lp.z,b3z = 3.-lp.z,b4z = 3.+lp.z;
	
	// Pre-compute expensive divisions
	double inv_size_x = 1. / (48.*lp.x);
	double inv_size_y = 1. / (48.*lp.y);
	double inv_size_z = 1. / (48.*lp.z);
	double oneTwelth = 1./12.;
	double inv_dx = 0;
	double inv_dy = 0;
	double inv_dz = 0;
	if (getDeriv) {
		inv_dx = 2.0 / GetDeltaX();
		inv_dy = 2.0 / GetDeltaY();
		inv_dz = 2.0 / GetDeltaZ();
	}
	
	int i=0;
	for(int id=0;id<64;id++)
	{	// x direction
		xp = fabs(xi->x - g3xii[id]);
		if(xp>=b4x) continue;
		
		// y direction
		yp = fabs(xi->y - g3eti[id]);
		if(yp>=b4y) continue;
		
		// y direction
		zp = fabs(xi->z - g3zti[id]);
		if(zp>=b4z) continue;
		
		// shape function i
		
		// x direction
		if(xp < b1x)
			Svpx = (9.-lp.x*lp.x-3*xp*xp)*oneTwelth;
		else if(xp < b2x)
		{	arg = xp-1.;
			double lp2 = lp.x*lp.x;
			double lp3 = lp2*lp.x;
			Svpx = (9.*lp2*arg + 3.*arg*arg*arg + 3.*lp.x*(15.-xp*(6.+xp)) - lp3)*inv_size_x;
		}
		else if (xp <= b3x)
		{	arg = xp-3.;
			Svpx = (lp.x*lp.x+3.*arg*arg)*0.5*oneTwelth;
		}
		else
		{	arg = 3. + lp.x - xp;
			Svpx = arg*arg*arg*inv_size_x;
		}
		
		// y direction
		if(yp < b1y)
			Svpy = (9.-lp.y*lp.y-3*yp*yp)*oneTwelth;
		else if(yp < b2y)
		{	arg = yp-1.;
			double lp2 = lp.y*lp.y;
			double lp3 = lp2*lp.y;
			Svpy = (9.*lp2*arg + 3.*arg*arg*arg + 3.*lp.y*(15.-yp*(6.+yp)) - lp3)*inv_size_y;
		}
		else if (yp <= b3y)
		{	arg = yp-3.;
			Svpy = (lp.y*lp.y+3.*arg*arg)*0.5*oneTwelth;
		}
		else
		{	arg = 3. + lp.y - yp;
			Svpy = arg*arg*arg*inv_size_y;
		}
		
		// z direction
		if(zp < b1z)
			Svpz = (9.-lp.z*lp.z-3*zp*zp)*oneTwelth;
		else if(zp < b2z)
		{	arg = zp-1.;
			double lp2 = lp.z*lp.z;
			double lp3 = lp2*lp.z;
			Svpz = (9.*lp2*arg + 3.*arg*arg*arg + 3.*lp.z*(15.-zp*(6.+zp)) - lp3)*inv_size_z;
		}
		else if (zp <= b3z)
		{	arg = zp-3.;
			Svpz = (lp.z*lp.z+3.*arg*arg)*0.5*oneTwelth;
		}
		else
		{	arg = 3. + lp.z - zp;
			Svpz = arg*arg*arg*inv_size_z;
		}

		// combine
		sfxn[i] = Svpx*Svpy*Svpz;
		
		// find shape function gradient at (xp,yp)
		
		if(getDeriv)
		{	xsign = xi->x>g3xii[id] ? 1. : -1.;
			ysign = xi->y>g3eti[id] ? 1. : -1.;
			zsign = xi->y>g3zti[id] ? 1. : -1.;
			
			// x gradient
			if(xp < b1x)
				dSvpx = -0.5*xp;
			else if(xp < b2x)
			{	arg = xp-1.;
				dSvpx = (3*lp.x*lp.x + 3.*arg*arg - 2.*lp.x*(3.+xp))*3.*inv_size_x;
			}
			else if (xp <= b3x)
				dSvpx = 0.25*(xp-3.);
			else
			{	arg = 3. + lp.x - xp;
				dSvpx = -arg*arg*3.*inv_size_x;
			}
			
			// y gradient
			if(yp < b1y)
				dSvpy = -0.5*yp;
			else if(yp < b2y)
			{	arg = yp-1.;
				dSvpy = (3*lp.y*lp.y + 3.*arg*arg - 2.*lp.y*(3.+yp))*3.*inv_size_y;
			}
			else if (yp <= b3y)
				dSvpy = 0.25*(yp-3.);
			else
			{	arg = 3. + lp.y - yp;
				dSvpy = -arg*arg*3.*inv_size_y;
			}
			
			// z gradient
			if(zp < b1z)
				dSvpz = -0.5*zp;
			else if(zp < b2z)
			{	arg = zp-1.;
				dSvpz = (3*lp.z*lp.z + 3.*arg*arg - 2.*lp.z*(3.+zp))*3.*inv_size_z;
			}
			else if (zp <= b3z)
				dSvpz = 0.25*(zp-3.);
			else
			{	arg = 3. + lp.z - zp;
				dSvpz = -arg*arg*3.*inv_size_z;
			}
			
			// combine
			xDeriv[i] = xsign*dSvpx*Svpy*Svpz*inv_dx;
			yDeriv[i] = ysign*Svpx*dSvpy*Svpz*inv_dy;
			zDeriv[i] = zsign*Svpx*Svpy*dSvpz*inv_dz;
		}
		
		// the node
		i++;
		nds[i] = nodes[0] + x3off[id]*mpmgrid.xplane + y3off[id]*mpmgrid.yplane + z3off[id]*mpmgrid.zplane;
	}
	
	// number of nodes found - may be has high as 16 and that is OK
	nds[0] = i;
}

#pragma mark EightNodeIsoparamBrick::Accessors

// element name as an ID
short EightNodeIsoparamBrick::ElementName(void) { return(EIGHT_NODE_ISO_BRICK); }

// number of nodes in this element
int EightNodeIsoparamBrick::NumberNodes(void) const { return 8; }

// Get x-y area, z thickness, or volumne - all orthogonal brick
double EightNodeIsoparamBrick::GetArea(void) const { return (xmax-xmin)*(ymax-ymin); }
double EightNodeIsoparamBrick::GetVolume(void) const { return (xmax-xmin)*(ymax-ymin)*(zmax-zmin); }
double EightNodeIsoparamBrick::GetThickness(void) const { return zmax-zmin; }

