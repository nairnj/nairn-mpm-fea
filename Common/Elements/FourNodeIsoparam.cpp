/******************************************************************************** 
    FourNodeIsoparam.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 24 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Elements/FourNodeIsoparam.hpp"
#include "Nodes/NodalPoint.hpp"
#ifdef MPM_CODE
	#include "NairnMPM_Class/MeshInfo.hpp"
    #include "NairnMPM_Class/NairnMPM.hpp"
	#include "MPM_Classes/MPMBase.hpp"
	#include <map>
	#include <vector>
#endif

// globals for node locations
static double xii[4]={-1.,1.,1.,-1.};
static double eti[4]={-1.,-1.,1.,1.};

// globals for GIMP node locations
#ifdef MPM_CODE
static double gxii[16]={-1.,1.,1.,-1.,-3.,-1.,1.,3.,3.,3.,3.,1.,-1.,-3.,-3.,-3.};
static double geti[16]={-1.,-1.,1.,1.,-3.,-3.,-3.,-3.,-1.,1.,3.,3.,3.,3.,1.,-1.};
// = int(0.5*(gxii[i]+1)) and int(0.5*(getai[i]+1))
static int xoff[16]={0,1,1,0,-1,0,1,2,2,2,2,1,0,-1,-1,-1};
static int yoff[16]={0,0,1,1,-1,-1,-1,-1,0,1,2,2,2,2,1,0};
#endif

#pragma mark FourNodeIsoparam::Constructors and Destructor

#ifdef MPM_CODE
/* Main MPM constructor - passes on to Linear 2D
	But also sets extra node to match starting node
*/
FourNodeIsoparam::FourNodeIsoparam(int eNum,int *eNode) : Linear2D(eNum,eNode)
{
    nodes[4]=eNode[0];			// needed in GetArea(), PtInElement()
	thickness=1.;
}

#else
/* Main FEA constructor - passes on to Linear 2D
	But also sets extra node to match starting node
*/
FourNodeIsoparam::FourNodeIsoparam(int eNum,int *eNode,int eMat,double eAng,double eThick) : 
            Linear2D(eNum,eNode,eMat,eAng,eThick)
{
    nodes[4]=eNode[0];			// needed in GetArea(), PtInElement()
}
#endif

#pragma mark FourNodeIsoparam::Methods

/*  Load shape functions and optionally their derivatives into arrays
    All arrays are 0 based (pass &(*)[1] to get 1-based results)
        
    if getDeriv==FALSE
        xiDeriv[], etaDeriv[], and eNodes[] not used (can be NULL)
    else if not BMATRIX
        etaDeriv[] and xiDeriv[] are the derivatives vs eta and xi
    else if BMATRIX
        etaDeriv[] and xiDeriv[] are derivatives vs x and y
        sfx[] will be Ni/2
        needs eNodes[] to be input
*/
void FourNodeIsoparam::ShapeFunction(Vector *xi,int getDeriv,
		double *sfxn,double *xiDeriv,double *etaDeriv,Vector *eNodes,
                double *outDetjac,double *outAsr,double *asbe) const
{
    double temp1,temp2,jac[3][3],detjac,asr;
    int i;
    
    // shape function or derivatives
    for(i=0;i<4;i++)
    {	temp1=(1.+xii[i]*xi->x)/4.;
        temp2=(1.+eti[i]*xi->y)/4.;
        sfxn[i]=4.*temp1*temp2;
        if(getDeriv)
        {	xiDeriv[i] = xii[i]*temp2;
            etaDeriv[i] = eti[i]*temp1;
        }
    }
    
    // Get B matrix elements for gradients
    if(getDeriv==BMATRIX)
    {	// Find Jacobian Matrix
        jac[1][1]=0.;
        jac[1][2]=0.;
        jac[2][1]=0.;
        jac[2][2]=0.;
        for(i=0;i<4;i++)
        {   jac[1][1]+=xiDeriv[i]*eNodes[i].x;
            jac[1][2]+=xiDeriv[i]*eNodes[i].y;
            jac[2][1]+=etaDeriv[i]*eNodes[i].x;
            jac[2][2]+=etaDeriv[i]*eNodes[i].y;
        }
        detjac=jac[1][1]*jac[2][2]-jac[1][2]*jac[2][1];
        temp1=jac[1][1]/detjac;
        jac[1][1]=jac[2][2]/detjac;
        jac[2][2]=temp1;
        jac[1][2]=-jac[1][2]/detjac;
        jac[2][1]=-jac[2][1]/detjac;
            
        // For radial position for axisymmetric analyses
        asr=0.;
        for(i=0;i<4;i++)
            asr=asr+sfxn[i]*eNodes[i].x;
                    
        /* Load J(-1)(1,1) ¶Ni/¶xi + J(-1)(1,2) ¶Ni/¶eta into xiDeriv[i]
                J(-1)(2,1) ¶Ni/¶xi + J(-1)(2,2) ¶Ni/¶eta into etaDeriv[i]
                Ni/r into asbe[i] */
        for(i=0;i<4;i++)
        {   temp1=xiDeriv[i];
            temp2=etaDeriv[i];
            xiDeriv[i]=jac[1][1]*temp1 + jac[1][2]*temp2;
            etaDeriv[i]=jac[2][1]*temp1 + jac[2][2]*temp2;
            if(asr!=0.)
                asbe[i]=sfxn[i]/asr;
            else
                asbe[i]=0.;
        }
        
        // return results
        if(outDetjac!=NULL) *outDetjac=detjac;
        if(outAsr!=NULL) *outAsr=asr;
    }
}

#ifdef FEA_CODE

// Take stress at four gauss points and map them to 4 nodes in this element
// by using coordinate system on the gauss points (numbered ccw as 1,3,4,2)
// as mapping to -1 to 1 coordinates.
// See FEA notes on stress extrapolation
// sgp[i][j] is stress j (1 to 4) at Gauss point i (1 to numGauss)
// se[i][j] is output stress j (1 to 4) at node i (1 to numnds) (externed variable)
void FourNodeIsoparam::ExtrapolateGaussStressToNodes(double sgp[][5])
{
	double gpt = 0.577350269189626;
	double at=1.+1./gpt;
	double bt=1.-1./gpt;
	double at2=at*at/4.;
	double bt2=bt*bt/4.;
	double atbt=at*bt/4.;
	
	double qe[5][5];
	qe[1][1]=at2;
	qe[1][2]=atbt;
	qe[1][3]=atbt;
	qe[1][4]=bt2;
	qe[2][1]=atbt;
	qe[2][2]=bt2;
	qe[2][3]=at2;
	qe[2][4]=atbt;
	qe[3][1]=bt2;
	qe[3][2]=atbt;
	qe[3][3]=atbt;
	qe[3][4]=at2;
	qe[4][1]=atbt;
	qe[4][2]=at2;
	qe[4][3]=bt2;
	qe[4][4]=atbt;
	
	// hard coded to 4 nodes and 4 Gauss points
	double temp;
	int i,j,k;
	for(i=1;i<=4;i++)
	{	for(j=1;j<=4;j++)
		{   temp=0.;
			for(k=1;k<=4;k++)
				temp+=qe[i][k]*sgp[k][j];
			se[i][j]=temp;
		}
	}
}

#endif

#ifdef MPM_CODE

// get shape functions and optionally derivatives wrt x and y, but derivatives only work
// if it is a regular array. Shape functions are general
// For axisymmetric MPM, make sure zDeriv is not NULL and load with shape function/rp
void FourNodeIsoparam::ShapeFunction(Vector *xi,int getDeriv,double *sfxn,
										double *xDeriv,double *yDeriv,double *zDeriv) const
{
    double temp1,temp2,dx,dy,rp;
	int i;
    
    if(getDeriv)
    {   dx = GetDeltaX();
        dy = GetDeltaY();
        rp = GetCenterX() + 0.5*xi->x*GetDeltaX();
    }
    
    // just shape function
    for(i=0;i<4;i++)
    {	temp1=(1.+xii[i]*xi->x);
        temp2=(1.+eti[i]*xi->y);
        sfxn[i]=0.25*temp1*temp2;
		if(getDeriv)
		{	xDeriv[i]=0.5*xii[i]*temp2/dx;
			yDeriv[i]=0.5*eti[i]*temp1/dy;
            if(fmobj->IsAxisymmetric()) zDeriv[i] = sfxn[i]/rp;
		}
    }
}

// get B2SPLINE shape functions and optionally derivatives wrt x and y, but derivatives only work
// if it is a regular array. Shape functions are general
// For axisymmetric MPM, make sure zDeriv is not NULL and load with shape function/rp
void FourNodeIsoparam::SplineShapeFunction(int *nds,Vector *xi,int getDeriv,double *sfxn,
										double *xDeriv,double *yDeriv,double *zDeriv) const
{
	// constants
	double inv_dx = 0.,inv_dy=0.,inv_rp=0.;
	if(getDeriv)
	{   inv_dx = 1./GetDeltaX();
		inv_dy = 1./GetDeltaY();
		inv_rp = 1./(GetCenterX() + 0.5*xi->x*GetDeltaX());
	}

	// get cell shape functions
	double sx,sy,arg,etai,netai,temp1,temp2;
	int i=0;
	for(int id=0;id<16;id++)
	{	// xi direction
		etai = xi->x - gxii[id];
		temp1=fabs(etai);
		if(temp1>=3.) continue;
		
		// y direction
		netai = xi->y - geti[id];
		temp2=fabs(netai);
		if(temp2>=3.) continue;
		
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
		
		// combine
		sfxn[i] = sx*sy;
		
		// Note that derivative is (2/dx)dN/detai and (2/dy)dN/dnetai
		if(getDeriv)
		{	// x direction
			if(temp1<=1.0)
				xDeriv[i] = -etai*sy*inv_dx;
			else if(etai>=0.)
				xDeriv[i] = 0.5*(etai-3)*sy*inv_dx;
			else
				xDeriv[i] = 0.5*(etai+3)*sy*inv_dx;
			
			// y direction
			if(temp2<=1.0)
				yDeriv[i] = -netai*sx*inv_dy;
			else if(netai>=0.)
				yDeriv[i] = 0.5*(netai-3)*sx*inv_dy;
			else
				yDeriv[i] = 0.5*(netai+3)*sx*inv_dy;
			
			// for axisymmetric
			if(fmobj->IsAxisymmetric()) zDeriv[i] = sfxn[i]*inv_rp;
		}
		
		// assign node number
		i++;
		nds[i] = nodes[0] + xoff[id]*mpmgrid.xplane + yoff[id]*mpmgrid.yplane;
	}
	
	// set number
	// this better be <=9 or bad things happen
	nds[0] = i;
}

//	Find extent of this element - called once at start (and must be called)
void FourNodeIsoparam::FindExtent(void)
{
	ElementBase::FindExtent();
	
	/* for speed in GetXiPos() calculations - if it a parallelogram
		Note: assumes element does not move. If it does, must recalculate these terms
	*/
	pgElement=PGRAM_ELEMENT;
	
	// are edges parallel wrt x coordinate?
	double xdel=fabs(nd[nodes[2]]->x-nd[nodes[1]]->x);
    double ydel=fabs(nd[nodes[3]]->x-nd[nodes[0]]->x);
	if(!DbleEqual(xdel,ydel) && xdel>1.e-10 && ydel>1.e-10)
	{	pgElement=QUAD_ELEMENT;
		return;
	}
	
	// are edges parallel wrt y coordinate?
	xdel=fabs(nd[nodes[2]]->y-nd[nodes[1]]->y);
	ydel=fabs(nd[nodes[3]]->y-nd[nodes[0]]->y);
	if(!DbleEqual(xdel,ydel) && xdel>1.e-10 && ydel>1.e-10)
	{	pgElement=QUAD_ELEMENT;
		return;
	}
	
	// is it a rectangle with [0] on lower left
	if(DbleEqual(nd[nodes[0]]->x,nd[nodes[3]]->x) && DbleEqual(nd[nodes[0]]->y,nd[nodes[1]]->y)
		&& (nd[nodes[1]]->x>nd[nodes[0]]->x) && (nd[nodes[3]]->y>nd[nodes[0]]->y))
	{	pgElement = RECT_ELEMENT;
		return;
	}

	// precalculate useful terms
	pgTerm[0]=nd[nodes[0]]->x+nd[nodes[1]]->x+nd[nodes[2]]->x+nd[nodes[3]]->x;
	pgTerm[1]=nd[nodes[1]]->x+nd[nodes[2]]->x-nd[nodes[0]]->x-nd[nodes[3]]->x;
	pgTerm[2]=nd[nodes[2]]->x+nd[nodes[3]]->x-nd[nodes[0]]->x-nd[nodes[1]]->x;
	pgTerm[3]=nd[nodes[0]]->y+nd[nodes[1]]->y+nd[nodes[2]]->y+nd[nodes[3]]->y;
	pgTerm[4]=nd[nodes[1]]->y+nd[nodes[2]]->y-nd[nodes[0]]->y-nd[nodes[3]]->y;
	pgTerm[5]=nd[nodes[2]]->y+nd[nodes[3]]->y-nd[nodes[0]]->y-nd[nodes[1]]->y;
	double det=pgTerm[1]*pgTerm[5]-pgTerm[2]*pgTerm[4];
	pgTerm[1]/=det;
	pgTerm[2]/=det;
	pgTerm[4]/=det;
	pgTerm[5]/=det;
}

/* Find dimensionless coordinates analytically if possible
	methods - on input, xi and eta should be initial
	guess in case needed for numerical
	
	Note: analytical possible if parallelogram, and terms precalculated
		for speed in FindExtent() above.
	
	Only used in MPM
*/
void FourNodeIsoparam::GetXiPos(const Vector *pos,Vector *xipos) const
{
	if(pgElement==RECT_ELEMENT)
	{	xipos->x=(2.*pos->x-xmin-xmax)/GetDeltaX();
		xipos->y=(2.*pos->y-ymin-ymax)/GetDeltaY();
	}
	else if(pgElement==PGRAM_ELEMENT)
	{	// analytical solution for parallelograms
		double xdel=4.*pos->x-pgTerm[0];
		double ydel=4.*pos->y-pgTerm[3];
		xipos->x=(pgTerm[5]*xdel-pgTerm[2]*ydel);
		xipos->y=(pgTerm[1]*ydel-pgTerm[4]*xdel);
	}
	else
		ElementBase::GetXiPos(pos,xipos);
}

// Find Cartesion position from natural coordinates
void FourNodeIsoparam::GetPosition(Vector *xipos,Vector *pos)
{
	if(pgElement==RECT_ELEMENT)
	{	pos->x = 0.5*(xmin+xmax+GetDeltaX()*xipos->x);
		pos->y = 0.5*(ymin+ymax+GetDeltaY()*xipos->y);
	}
	else if(pgElement==PGRAM_ELEMENT)
	{	double xarg = 0.25*(xipos->x + pgTerm[5]*pgTerm[0] - pgTerm[2]*pgTerm[3]);
		double yarg = 0.25*(xipos->y + pgTerm[1]*pgTerm[3] - pgTerm[4]*pgTerm[0]);
		double denom = (pgTerm[1]*pgTerm[5] - pgTerm[2]*pgTerm[4]);
		pos->x = (pgTerm[1]*xarg + pgTerm[2]*yarg) / denom;
		pos->y = (pgTerm[4]*xarg + pgTerm[5]*yarg) / denom;
		pos->z = 0.;
	}
	else
		ElementBase::GetPosition(xipos,pos);
}

// see if this element is rectangle in cartesion coordinates returning TRUE or FALSE
// if true, dx and dy set to element dimensions and dz set to zero
int FourNodeIsoparam::Orthogonal(double *dx,double *dy,double *dz)
{
	int i;
	double xdel,ydel;
	
	*dx=0.;
	*dy=0.;
	*dz=0.;
	for(i=0;i<3;i++)
    {	xdel=fabs(nd[nodes[i+1]]->x-nd[nodes[i]]->x);
		ydel=fabs(nd[nodes[i+1]]->y-nd[nodes[i]]->y);
		if(xdel>=1.e-12 && ydel>=1.e-12) return FALSE;
		*dx=fmax(xdel,*dx);
		*dy=fmax(ydel,*dy);
	}
    xdel=fabs(nd[nodes[0]]->x-nd[nodes[3]]->x);
	ydel=fabs(nd[nodes[0]]->y-nd[nodes[3]]->y);
	if(xdel>=1.e-12 && ydel>=1.e-12) return FALSE;
	*dx=fmax(xdel,*dx);
	*dy=fmax(ydel,*dy);
	return TRUE;
}

// Faster if rectangular element
short FourNodeIsoparam::PtInElement(Vector &pt) const
{
	if(pgElement==RECT_ELEMENT)
	{	if(pt.x<xmin || pt.x>=xmax) return false;
		if(pt.y<ymin || pt.y>=ymax) return false;
		return true;
	}
	
	// pass to subclass
	return Linear2D::PtInElement(pt);
}

#pragma mark FourNodeIsoparam::GIMP Methods

// get GIMP shape functions and optionally derivatives wrt x and y
// assumed to be properly numbered regular array
// input *xi position in element coordinate
// output number of nodes in nds[0] and nodes in nds[1] to nds[nds[0]
void FourNodeIsoparam::GimpShapeFunction(Vector *xi,int *nds,int getDeriv,double *sfxn,
										 double *xDeriv,double *yDeriv,double *zDeriv,Vector &lp) const
{
	int i;
	double xp,yp,Svpx,Svpy,dSvpx,dSvpy,xsign,ysign,argx=0.,argy=0.;
	
	// L is the cell spacing, 2*lpi is the current particle size (dimensionless range -1 to 1).
	// The deformation of the particle is not considered yet.
    double q1x = 2.-lp.x, q2x = 2.+lp.x;
    double q1y = 2.-lp.y, q2y = 2.+lp.y;
	
	// Pre-compute expensive divisions
	double inv_size_x = 1. / (4.*lp.x);
	double inv_size_y = 1. / (4.*lp.y);
	double inv_dx = 0; 
	double inv_dy = 0;
	if (getDeriv) {
		 inv_dx = 2.0 / GetDeltaX();
		 inv_dy = 2.0 / GetDeltaY();
	}
	
	i=0;
	for(int id=0;id<16;id++)
	{	// x direction
		xp = fabs(xi->x - gxii[id]);
		if(xp>=q2x) continue;
		
		// y direction
		yp = fabs(xi->y - geti[id]);
		if(yp>=q2y) continue;
		
		// shape function for node i
		
		// x direction
		if (xp < lp.x)
			Svpx = ((4. - lp.x)*lp.x - xp*xp)*inv_size_x; 		// if lp=0.5: -(4.*xp*xp-7.)/8.;
		else if(xp<=q1x)
			Svpx = (2.-xp)/2.;
		else
		{	argx = (q2x - xp)*inv_size_x; 						// if lp=0.5: (5.-2.*xp)/4
			Svpx = 2.*lp.x*argx*argx;							// if lp=0.5: argx*argx
		}
		
		// y direction
		if (yp < lp.y)
			Svpy = ((4. - lp.y)*lp.y - yp*yp)*inv_size_y; 		// if lp=0.5: -(4.*yp*yp-7.)/8.;
		else if(yp<=q1y)
			Svpy = (2.-yp)/2.;
		else
		{	argy = (q2y - yp)*inv_size_y; 						// if lp=0.5: (5.-2.*yp)/4
			Svpy = 2.*lp.y*argy*argy;							// if lp=0.5: (5.-2.*yp)^2/16
		}
		
		sfxn[i] = Svpx*Svpy;
		
		// find shape function gradient at (xp,yp)
		
		if(getDeriv)
		{	xsign = xi->x>gxii[id] ? 1. : -1.;
			ysign = xi->y>geti[id] ? 1. : -1.;
			
			// x gradient
			if (xp < lp.x)
				dSvpx = -xp*inv_size_x*2.0; 					// if lp=0.5: -xp
			else if(xp<=q1x)
				dSvpx = -0.5;
			else
				dSvpx = -argx;
			
			// y gradient
			if (yp < lp.y)
				dSvpy = -yp*inv_size_y*2.0; 					// if lp=0.5: -xp
			else if(yp<=q1y)
				dSvpy = -0.5;
			else
				dSvpy = -argy;
			
			xDeriv[i] = xsign*dSvpx*Svpy*inv_dx;
			yDeriv[i] = ysign*Svpx*dSvpy*inv_dy;
		}
		
		// assign node
		i++;
		nds[i] = nodes[0] + xoff[id]*mpmgrid.xplane + yoff[id]*mpmgrid.yplane;
	}
	
	// number of nodes (better be less than space available in nds[])
	nds[0] = i;
}

// get GIMP shape functions and optionally derivatives wrt x and y
// assumed to be properly numbered regular array
// input *xi position in element coordinate
// output number of nodes in nds[0] and nodes in nds[1] to nds[nds[0]]
void FourNodeIsoparam::GimpShapeFunctionAS(Vector *xi,int *nds,int getDeriv,double *sfxn,
										   double *xDeriv,double *yDeriv,double *zDeriv,Vector &lp) const
{
	int i,n;
	double xp,yp,ri,nr,Svpx,Svpy,dSvpx,dSvpy,pTr,ysign,argx=0.,argy=0.;
	
	// L is the cell spacing, 2*lp is the current particle size.
	// assuming the particle size is the same in x and y direction in element coordinates
	// the deformation of the particle is not considered here.
	// assumes 4 particles per element
	double q2x = 2.-lp.x, q3x = 2.+lp.x;
	double q2y = 2.-lp.y, q3y = 2.+lp.y;
	double dx = GetDeltaX();
	double midx = GetCenterX();
	double dy = GetDeltaY();

	// Pre-compute expensive divisions
	double inv_size_x = 1. / (4.*lp.x);
	double inv_size_y = 1. / (4.*lp.y);
	double inv_dx = 0;
	double inv_dy = 0;
	if (getDeriv) {
		inv_dx = 2.0 / dx;
		inv_dy = 2.0 / dy;
	}

	i = 0;
	for(int id=0;id<16;id++)
	{	// x direction
		xp = xi->x-gxii[id];			// signed xp
		if(fabs(xp)>=q3x) continue;
		
		// y direction
		yp = fabs(xi->y-geti[id]);		// first quadrant for yp>0
		if(yp>=q3y) continue;
		
		// find nodal position based on node numbers and nodal column number
		ri = midx+0.5*gxii[id]*dx;		// -1.5dx, -.5dx, +.5dx, 1.5dx  from element center
		nr = ri/dx;						// radial node number (should be an integer)
		
		// truncate to zero if node at negative r value
		if(nr<-0.5) continue;
		
		// shape function i
		
		// truncate near r=0
		if(nr<0.5)
		{	// node at r=0 (note that xp<0 never occurs and no need to check)
			n=0;
		}
		else if(nr<1.5)
		{	// node at r=dx *note that xp<-2 never occurs and no need to check)
			n=1;
		}
		else
		{	// node at r>=2dx
			n=2;
		}
		
		// X direction with truncation
		// note argx used below and only set for xp<-q2x (n>=2) and xp>q2x
		if(xp<-q2x)
		{	// note that n==0 never occurs here
			if(n==1)
				Svpx = (q3x+xp)/3.;						// note that xp>-2 always (no need to check)
			else
			{	argx = (q3x+xp)*inv_size_x;
				Svpx = 2.*lp.x*argx*argx*(1.-(2.*(1.-lp.x)+xp)/(3.*(2.*nr+xp)));
			}
		}
		else if(xp<=-lp.x)
			Svpx = 0.5*(2.+xp) + lp.x*lp.x/(6.*(2.*nr+xp));
		else if(xp<lp.x)
		{	if(n==0)
				Svpx = (3.-lp.x-xp)/3.;					// note that xp>0 always (no need to check)
			else
				Svpx = ((4.-lp.x)*lp.x-xp*xp)*inv_size_x + xp*(xp*xp-3.*lp.x*lp.x)/(12.*lp.x*(2.*nr+xp));
		}
		else if(xp<=q2x)
			Svpx = 0.5*(2.-xp) - lp.x*lp.x/(6.*(2.*nr+xp));		// if lp=0.5: (2.-xp)/2. - 1./(24.*(2.*nr+xp));
		else
		{	argx = (q3x-xp)*inv_size_x;
			Svpx = 2.*lp.x*argx*argx*(1.+(2.*(1.-lp.x)-xp)/(3.*(2.*nr+xp)));
		}
		
		// Y direction like planar GIMP
		// note argy used below and only set for yp>q2y
		if(yp<lp.y)
			Svpy = ((4.-lp.y)*lp.y-yp*yp)*inv_size_y;
		else if(yp<=q2y)
			Svpy = 0.5*(2.-yp);
		else
		{	argy = (q3y-yp)*inv_size_y;
			Svpy = 2.*lp.y*argy*argy;
		}
		
		sfxn[i] = Svpx*Svpy;
		
		// find shape function at (xp,yp)
		if(getDeriv)
		{	ysign = xi->y>geti[id] ? 1. : -1.;
			
			// X part with truncation options
			if(xp<-q2x)
			{	// note that n==0 never occurs here
				if(n==1)
					dSvpx = 0.5;                    // note that xp>-2 always (no need to check)
				else
					dSvpx = argx*(1 - (q2x+xp)/(2.*(2.*nr+xp)));
			}
			else if(xp<=-lp.x)
				dSvpx = 0.5;
			else if(xp<lp.x)
			{	if(n==0)
					dSvpx = -0.5;
				else
					dSvpx = -2.0*xp*inv_size_x - (lp.x*lp.x-xp*xp)/(4.*lp.x*(2.*nr+xp));
			}
			else if(xp<=q2x)
				dSvpx = -0.5;                    // note that xp>0 always (no need to check)
			else
				dSvpx = -argx*(1. + (q2x-xp)/(2.*(2.*nr+xp)));
			
			// Y part
			if(yp<lp.y)
				dSvpy = -ysign*yp/(2.*lp.y);
			else if(yp<=q2y)
				dSvpy = -0.5*ysign;
			else
				dSvpy = -argy*ysign;
			
			xDeriv[i] = dSvpx*Svpy*inv_dx;
			yDeriv[i] = Svpx*dSvpy*inv_dy;
			
			// Find plane Sx/(2n+xp)
			if(xp<-q2x)
			{	// note that n==0 never occurs here
				if(n==1)
					pTr = 0.5;
				else
				{	argx = (q3x+xp)*inv_size_x;
					pTr = 2.*lp.x*argx*argx/(2.*nr+xp);
				}
			}
			else if(xp<=-lp.x)
				pTr = (2.+xp)/(2.*(2.*nr+xp));              // which is 0.5 for n==1
			else if(xp<lp.x)
			{	if(n==0)
					pTr = 2./(xp+lp.x) - 0.5;
				else
					pTr = ((4.-lp.x)*lp.x-xp*xp)/(4.*lp.x*(2.*nr+xp));
			}
			else if(xp<=q2x)
				pTr = (2.-xp)/(2.*(2.*nr+xp));
			else
			{	argx = (q3x-xp)*inv_size_x;
				pTr = 2.*lp.x*argx*argx/(2.*nr+xp);
			}
			
			zDeriv[i] = pTr*Svpy*inv_dx;
		}
		
		// get node
		i++;
		nds[i] = nodes[0] + xoff[id]*mpmgrid.xplane + yoff[id]*mpmgrid.yplane;
	}
	
	// total number (better be <=9)
	nds[0] = i;
}

// get B2GIMP shape functions and optionally derivatives wrt x and y
// assumed to be properly numbered regular array
// input *xi position in element coordinate and ndIDs[0]... is which nodes (0-15)
void FourNodeIsoparam::BGimpShapeFunction(Vector *xi,int *nds,int getDeriv,double *sfxn,
										 double *xDeriv,double *yDeriv,double *zDeriv,Vector &lp) const
{
	double xp,yp,Svpx,Svpy,dSvpx,dSvpy,xsign,ysign,arg;
	
	// L is the cell spacing, 2*lpi is the current particle size (dimensionless range -1 to 1).
	// Breakpoints on positive side of the node
	double b1x = 1.-lp.x,b2x = 1.+lp.x,b3x = 3.-lp.x,b4x = 3.+lp.x;
	double b1y = 1.-lp.y,b2y = 1.+lp.y,b3y = 3.-lp.y,b4y = 3.+lp.y;
	
	// Pre-compute expensive divisions
	double inv_size_x = 1. / (48.*lp.x);
	double inv_size_y = 1. / (48.*lp.y);
	double oneTwelth = 1./12.;
	double inv_dx = 0;
	double inv_dy = 0;
	if (getDeriv) {
		inv_dx = 2.0 / GetDeltaX();
		inv_dy = 2.0 / GetDeltaY();
	}
	double lpx2=lp.x*lp.x;
	double lpy2=lp.y*lp.y;

	int i=0;
	for(int id=0;id<16;id++)
	{	// x direction
		xp = fabs(xi->x - gxii[id]);
		if(xp>=b4x) continue;
		
		// y direction
		yp = fabs(xi->y - geti[id]);
		if(yp>=b4y) continue;
		
		// shape function i
		
		if(xp < b1x)
			Svpx = (9.-lpx2-3.*xp*xp)*oneTwelth;
		else if(xp < b2x)
		{	arg = xp-1.;
			Svpx = (lpx2*(9.*arg-lp.x) + 3.*arg*arg*arg + 3.*lp.x*(15.-xp*(6.+xp)))*inv_size_x;
		}
		else if (xp <= b3x)
		{	arg = xp-3.;
			Svpx = (lpx2+3.*arg*arg)*0.5*oneTwelth;
		}
		else
		{	arg = 3. + lp.x - xp;
			Svpx = arg*arg*arg*inv_size_x;
		}
		
		if(yp < b1y)
			Svpy = (9.-lpy2-3.*yp*yp)*oneTwelth;
		else if(yp < b2y)
		{	arg = yp-1.;
			Svpy = (lpy2*(9.*arg-lp.y) + 3.*arg*arg*arg + 3.*lp.y*(15.-yp*(6.+yp)))*inv_size_y;
		}
		else if (yp <= b3y)
		{	arg = yp-3.;
			Svpy = (lpy2+3.*arg*arg)*0.5*oneTwelth;
		}
		else
		{	arg = 3. + lp.y - yp;
			Svpy = arg*arg*arg*inv_size_y;
		}
		
		sfxn[i] = Svpx*Svpy;
		
		// find shape function gradient at (xp,yp)
		
		if(getDeriv)
		{	xsign = xi->x>gxii[id] ? 1. : -1.;
			ysign = xi->y>geti[id] ? 1. : -1.;
			
			if(xp < b1x)
				dSvpx = -0.5*xp;
			else if(xp < b2x)
			{	arg = xp-1.;
				dSvpx = (3*lpx2 + 3.*arg*arg - 2.*lp.x*(3.+xp))*3.*inv_size_x;
			}
			else if (xp <= b3x)
				dSvpx = 0.25*(xp-3.);
			else
			{	arg = 3. + lp.x - xp;
				dSvpx = -arg*arg*3.*inv_size_x;
			}
			
			if(yp < b1y)
				dSvpy = -0.5*yp;
			else if(yp < b2y)
			{	arg = yp-1.;
				dSvpy = (3*lpy2 + 3.*arg*arg - 2.*lp.y*(3.+yp))*3.*inv_size_y;
			}
			else if (yp <= b3y)
				dSvpy = 0.25*(yp-3.);
			else
			{	arg = 3. + lp.y - yp;
				dSvpy = -arg*arg*3.*inv_size_y;
			}
			
			xDeriv[i] = xsign*dSvpx*Svpy*inv_dx;
			yDeriv[i] = ysign*Svpx*dSvpy*inv_dy;
		}
		
		// the node
		i++;
		nds[i] = nodes[0] + xoff[id]*mpmgrid.xplane + yoff[id]*mpmgrid.yplane;
	}
	
	// number of nodes found - may be has high as 16 and that is OK
	nds[0] = i;
}

// get B2GIMP shape functions and optionally derivatives wrt x and y
// assumed to be properly numbered regular array
// input *xi position in element coordinate and ndIDs[0]... is which nodes (0-15)
void FourNodeIsoparam::BGimpShapeFunctionAS(Vector *xi,int *nds,int getDeriv,double *sfxn,
										  double *xDeriv,double *yDeriv,double *zDeriv,Vector &lp) const
{
	int i;
	double xp,yp,Svpx,Svpy,dSvpx,dSvpy,ysign,arg,ri,nr,Tr;
	
	// L is the cell spacing, 2*lpi is the current particle size (dimensionless range -1 to 1).
	// Breakpoints on positive side of the node
	double b1x = 1.-lp.x,b2x = 1.+lp.x,b3x = 3.-lp.x,b4x = 3.+lp.x;
	double b1y = 1.-lp.y,b2y = 1.+lp.y,b3y = 3.-lp.y,b4y = 3.+lp.y;
	
	// Pre-compute expensive divisions
	double inv_size_x = 1. / (48.*lp.x);
	double inv_size_y = 1. / (48.*lp.y);
	double oneTwelth = 1./12.;
	double dx = GetDeltaX();
	double midx = GetCenterX();
	double dy = GetDeltaY();
	double inv_dx = 0;
	double inv_dy = 0;
	if (getDeriv) {
		inv_dx = 2.0 / dx;
		inv_dy = 2.0 / dy;
	}
	double lpx2 = lp.x*lp.x;
	double lpy2 = lp.y*lp.y;

	i = 0;
	for(int id=0;id<16;id++)
	{	// x direction
		xp = xi->x - gxii[id];
		if(fabs(xp)>=b4x) continue;
		
		// y direction
		yp = fabs(xi->y - geti[id]);
		if(yp>=b4y) continue;
		
		// find nodal position based on node numbers and nodal column number
		ri = midx+0.5*gxii[id]*dx;		// -1.5dx, -.5dx, +.5dx, 1.5dx  from element center
		nr = ri/dx;						// radial node number (should be an integer when r<=2*dx in grid)
		
		// shape function i
		// Particle radial position xp+2nr and radial extent xp-lp.x+2nr to xp+lp.x+2nr
		// Because xp+2nr>=0, then always have xp > -2nr
		// truncate if left edge<0 or xp < lp.x-2nr
		
		// Handle special case nodes first
		bool truncated = false;
		if(nr<-1.5)
		{	// node at r=-2dx should never occur
			continue;
		}
		else if(nr<-0.5)
		{	// nodes at r = -dx (n=-1), is left edge<0
			if(xp<2.+lp.x)
			{	arg = lp.x+xp-2.;
				if(xp<b3x)
					Svpx = 0.25*(6.-arg*(8.-3.*arg))*oneTwelth;
				else
					Svpx = 1./(48.*arg*arg);
				truncated = true;
			}
		}
		else if(nr<0.5)
		{	// nodes at r = 0 (n=0), is left edge<0
			if(xp<lp.x)
			{	arg = lp.x+xp;
				if(xp<=b1x)
					Svpx = 0.125*(6.-arg*arg);
				else
					Svpx = 0.0625*(18.-arg*(8.-arg)-1./(arg*arg));
				truncated = true;
			}
		}
		else if(nr<1.5)
		{	// nodes at r = dx (n=1), is left edge<0
			if(xp<-2.+lp.x)
			{	arg = lp.x+xp+2.;
				if(xp<=-b2x)
					Svpx = 0.25*(6.+arg*(8.+3.*arg))*oneTwelth;
				else
					Svpx = 2.*arg/3. - 0.0625*(4.+2.*arg*arg-1./(arg*arg));
				truncated = true;
			}
		}
		
		double twonix = 1./(2.*nr+xp);
		if(!truncated)
		{	// Seven pieces in function
			if(xp<-b3x)
			{	arg = 3. + lp.x + xp;
				Tr = arg*arg*arg*inv_size_x;
				Svpx = Tr*(1. - 0.25*(3.*(1.-lp.x)+xp)*twonix);
			}
			else if(xp<-b2x)
			{	arg = xp+3.;
				Tr = 0.5*(lpx2+3.*arg*arg)*oneTwelth;
				Svpx = Tr + lpx2*arg*twonix*oneTwelth;
			}
			else if(xp<-b1x)
			{	arg = xp+1.;
				double arg2 = arg*arg;
				Tr = (-lpx2*(9.*arg +lp.x) - 3.*arg2*arg + 3.*lp.x*(15.+xp*(6.-xp)))*inv_size_x;
				Svpx = Tr + 0.25*(3.*arg2*(arg2-6.*lpx2) - lpx2*lp.x*(9.*lp.x+8.*(xp-3.)))*twonix*inv_size_x;
			}
			else if(xp<b1x)
			{	Tr = (9.-lpx2 - 3.*xp*xp)*oneTwelth;
				Svpx = Tr - 2.*lpx2*xp*twonix*oneTwelth;
			}
			else if(xp<b2x)
			{	arg = xp-1.;
				double arg2 = arg*arg;
				Tr = (lpx2*(9.*arg-lp.x) + 3.*arg2*arg + 3.*lp.x*(15.-xp*(6.+xp)))*inv_size_x;
				Svpx = Tr + 0.25*(3.*arg2*(6.*lpx2-arg2) + lpx2*lp.x*(9.*lp.x-8.*(xp+3.)))*twonix;
			}
			else if(xp<b3x)
			{	arg = xp-3.;
				Tr = 0.5*(lpx2+3.*arg*arg)*oneTwelth;
				Svpx = Tr + lpx2*arg*twonix*oneTwelth;
			}
			else
			{	arg = 3. + lp.x - xp;
				Tr = arg*arg*arg*inv_size_x;
				Svpx = Tr*(1. + 0.25*(3.*(1.-lp.x)-xp)*twonix);
			}
		}

		if(yp < b1y)
			Svpy = (9.-lpy2-3.*yp*yp)*oneTwelth;
		else if(yp < b2y)
		{	arg = yp-1.;
			Svpy = (lpy2*(9.*arg-lp.y) + 3.*arg*arg*arg + 3.*lp.y*(15.-yp*(6.+yp)))*inv_size_y;
		}
		else if (yp <= b3y)
		{	arg = yp-3.;
			Svpy = (lpy2+3.*arg*arg)*0.5*oneTwelth;
		}
		else
		{	arg = 3. + lp.y - yp;
			Svpy = arg*arg*arg*inv_size_y;
		}
		
		sfxn[i] = Svpx*Svpy;
		
		if(getDeriv)
		{	// Handle truncated nodes first
			if(truncated)
			{	// if handle is true from above, then not need to check range again
				if(nr<-0.5)
				{	// nodes at r = -dx (n=-1)
					arg = lp.x+xp-2.;
					if(xp<b3x)
					{	dSvpx = (2.*arg - 3.)*oneTwelth;
						Tr = (3./arg-3.+arg)*oneTwelth;
					}
					else
					{	dSvpx = -oneTwelth/(arg*arg);
						Tr = -dSvpx;
					}
				}
				else if(nr<0.5)
				{	// nodes at r = 0 (n=0)
					arg=lp.x+xp;
					if(xp<=b1x)
					{	dSvpx = -arg/3.;
						Tr = (9./arg-arg)/6.;
					}
					else
					{	dSvpx = (3. - arg*arg*(9. - 2.*arg))/(12.*arg*arg);
						Tr = (-3.+arg*(27.-arg*(9.-arg)))/(12.*arg*arg);
					}
				}
				else
				{	// nodes at r = dx (n=1)
					arg=lp.x+xp+2.;
					if(xp<=-b2x)
					{	dSvpx = (2.*arg + 3.)*oneTwelth;
						Tr = (3./arg+3.+arg)*oneTwelth;
					}
					else
					{	dSvpx = -(3.-arg*arg*(12.-4.*arg))/(12.*arg*arg);
						Tr = (3.-arg*(6.-2.*arg*(6.-arg)))/(12.*arg*arg);
					}
				}
			}
			else
			{	if(xp<-b3x)
				{	arg = 3. + lp.x + xp;
					dSvpx = arg*arg*inv_size_x*(3.-(3.-2.*lp.x+xp)*twonix);
				}
				else if(xp<-b2x)
				{	dSvpx = 0.25*(xp+3.) + lpx2*oneTwelth*twonix;
				}
				else if(xp<-b1x)
				{	arg = xp+1.;
					double arg2=arg*arg;
					dSvpx = -(3.*(3*lpx2 + 3.*arg2 - 2.*lp.x*(3.-xp))
									+(lpx2*(2.*lp.x+9.*arg) - 3.*arg2*arg)*twonix)*inv_size_x;
				}
				else if(xp<b1x)
				{	dSvpx = -0.5*xp - 2.*lpx2*oneTwelth*twonix;
				}
				else if(xp<b2x)
				{	arg = xp-1.;
					double arg2 = arg*arg;
					dSvpx = (3.*(3*lpx2 + 3.*arg2 - 2.*lp.x*(3.+xp))
								- (lpx2*(2.*lp.x-9.*arg) + 3.*arg*arg2)*twonix)*inv_size_x;
				}
				else if(xp<b3x)
				{	dSvpx = 0.25*(xp-3.) + lpx2*oneTwelth*twonix;
				}
				else
				{	arg = 3. + lp.x - xp;
					dSvpx = -arg*arg*inv_size_x*(3.+(3.-2.*lp.x-xp)*twonix);
				}
				
				// scale Tr by 1/r
				Tr *= twonix;
			}
			
			if(yp < b1y)
				dSvpy = -0.5*yp;
			else if(yp < b2y)
			{	arg = yp-1.;
				dSvpy = (3*lpy2 + 3.*arg*arg - 2.*lp.y*(3.+yp))*3.*inv_size_y;
			}
			else if (yp <= b3y)
				dSvpy = 0.25*(yp-3.);
			else
			{	arg = 3. + lp.y - yp;
				dSvpy = -arg*arg*3.*inv_size_y;
			}
			
			xDeriv[i] = dSvpx*Svpy*inv_dx;
			ysign = xi->y>geti[id] ? 1. : -1.;
			yDeriv[i] = ysign*Svpx*dSvpy*inv_dy;
			zDeriv[i] = Tr*Svpy*inv_dx;
			if(Tr>2.) cout << Tr << "," << truncated << "," << nr << endl;
		}
		
		// the node
		i++;
		nds[i] = nodes[0] + xoff[id]*mpmgrid.xplane + yoff[id]*mpmgrid.yplane;
	}
	
	// number of nodes found - may be has high as 16 and that is OK
	nds[0] = i;
}

// Get shape functions for exact traction calculations in 2D MPM
void FourNodeIsoparam::GridTractionFunction(Vector *xi1,Vector *xi2,bool isAxisymmetric,double *tfxn,int *nds,double rmid,double dr) const
{
	// get the nodes (1 based)
	int numnds;
	GetNodes(&numnds,nds);
	nds[0] = numnds;
	
	// get shape functions (0 based)
	for(int i=0;i<4;i++)
	{	tfxn[i] = 1. + 0.5*xii[i]*(xi1->x+xi2->x) + 0.5*eti[i]*(xi1->y+xi2->y)
		+ (xii[i]*eti[i]/6.)*(2.*(xi1->x*xi1->y+xi2->x*xi2->y) + xi1->x*xi2->y + xi2->x*xi1->y);
		if(isAxisymmetric)
		{	tfxn[i] = rmid*tfxn[i] + dr*( xii[i]*(xi2->x-xi1->x) + eti[i]*(xi2->y-xi1->x)
										 + xii[i]*eti[i]*(xi2->x*xi2->y - xi1->x*xi1->y) );
		}
	}
}

#endif		// MPM_CODE

#pragma mark FourNodeIsoparam::Accessors

// element name as an ID
short FourNodeIsoparam::ElementName(void) { return(FOUR_NODE_ISO); }

// number of nodes in this element
int FourNodeIsoparam::NumberNodes(void) const { return 4; }
