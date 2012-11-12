/******************************************************************************** 
    FourNodeIsoparam.cpp
    NairnMPM
    
    Created by John Nairn on Wed Jan 24 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Elements/FourNodeIsoparam.hpp"
#include "Nodes/NodalPoint.hpp"
#ifdef MPM_CODE
	#include "NairnMPM_Class/MeshInfo.hpp"
    #include "NairnMPM_Class/NairnMPM.hpp"      // +AS
#endif

// globals for node locations
static double xii[4]={-1.,1.,1.,-1.};
static double eti[4]={-1.,-1.,1.,1.};

// globals for GIMP node locations
#ifdef MPM_CODE
static double gxii[16]={-1.,1.,1.,-1.,-3.,-1.,1.,3.,3.,3.,3.,1.,-1.,-3.,-3.,-3.};
static double geti[16]={-1.,-1.,1.,1.,-3.,-3.,-3.,-3.,-1.,1.,3.,3.,3.,3.,1.,-1.};
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
                double *outDetjac,double *outAsr,double *asbe)
{
    double temp1,temp2,jac[3][3],detjac,asr;
    int i;
    
    // shape function or derivatives
    for(i=0;i<4;i++)
    {	temp1=(1.+xii[i]*xi->x)/4.;
        temp2=(1.+eti[i]*xi->y)/4.;
        sfxn[i]=4.*temp1*temp2;
        if(getDeriv)
        {   xiDeriv[i]=xii[i]*temp2;
            etaDeriv[i]=eti[i]*temp1;
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
										double *xDeriv,double *yDeriv,double *zDeriv)
{
    double temp1,temp2,dx,dy,xp;
	int i;
    
    if(getDeriv)
    {   dx = GetDeltaX();
        dy = GetDeltaY();
        xp = GetCenterX() + 0.5*xi->x*GetDeltaX();
    }
    
    // just shape function
    for(i=0;i<4;i++)
    {	temp1=(1.+xii[i]*xi->x);
        temp2=(1.+eti[i]*xi->y);
        sfxn[i]=0.25*temp1*temp2;
		if(getDeriv)
		{	xDeriv[i]=0.5*xii[i]*temp2/dx;
			yDeriv[i]=0.5*eti[i]*temp1/dy;
            if(fmobj->IsAxisymmetric()) zDeriv[i] = sfxn[i]/xp;
		}
    }
}

//	Find extent of this element - called once at start (and must be called)
void FourNodeIsoparam::FindExtent(void)
{
	ElementBase::FindExtent();
	
	/* for speed in GetXiPos() calculations - if it a parallelogram
		Note: assumes element does not move. If it does, must recalculate these terms
	*/
	pgElement=TRUE;
	
	// are edges parallel wrt x coordinate?
	double xdel=fabs(nd[nodes[2]]->x-nd[nodes[1]]->x);
    double ydel=fabs(nd[nodes[3]]->x-nd[nodes[0]]->x);
	if(!DbleEqual(xdel,ydel) && xdel>1.e-10 && ydel>1.e-10) pgElement=FALSE;
	
	if(pgElement)
	{	// are edges parallel wrt y coordinate?
		xdel=fabs(nd[nodes[2]]->y-nd[nodes[1]]->y);
		ydel=fabs(nd[nodes[3]]->y-nd[nodes[0]]->y);
		if(!DbleEqual(xdel,ydel) && xdel>1.e-10 && ydel>1.e-10) pgElement=FALSE;
	}

	// precalculate useful terms
	if(pgElement)
	{	pgTerm[0]=nd[nodes[0]]->x+nd[nodes[1]]->x+nd[nodes[2]]->x+nd[nodes[3]]->x;
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
}

/* Find dimensionless coordinates analytically if possible
	methods - on input, xi and eta should be initial
	guess in case needed for numerical
	
	Note: analytical possible if parallelogram, and terms precalculated
		for speed in FindExtent() above.
	
	Only used in MPM
*/
void FourNodeIsoparam::GetXiPos(Vector *pos,Vector *xipos)
{
	// analytical solution for parallelograms
	if(pgElement)
	{	double xdel=4.*pos->x-pgTerm[0];
		double ydel=4.*pos->y-pgTerm[3];
		xipos->x=(pgTerm[5]*xdel-pgTerm[2]*ydel);
		xipos->y=(pgTerm[1]*ydel-pgTerm[4]*xdel);
	}
	else
		ElementBase::GetXiPos(pos,xipos);
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

#pragma mark FourNodeIsoparam::GIMP Methods

// Get GIMP nodes around an element #num, but only where shape functions is non zero
// assumed to be properly numbered regular array
// load nodes into nds[1]... and node IDs (0-15) into ndIDs[0]...
void FourNodeIsoparam::GetGimpNodes(int *numnds,int *nds,int *ndIDs,Vector *xipos)
{
	// quadrant barriers assuming 4 particles
	double q1=-0.5,q2=0.5;
	
	// nodes directly associated with the element
	nds[1]=nodes[0];
	nds[2]=nodes[1];
	nds[3]=nodes[2];
	nds[4]=nodes[3];
	ndIDs[0]=0;
	ndIDs[1]=1;
	ndIDs[2]=2;
	ndIDs[3]=3;
	
	// lower y quadrant
	if(xipos->y<q1)
	{	if(xipos->x<q1)
		{	nds[5]=nodes[0]-mpmgrid.xplane-mpmgrid.yplane;
			nds[6]=nds[5]+mpmgrid.xplane;
			nds[7]=nds[6]+mpmgrid.xplane;
			nds[8]=nodes[3]-mpmgrid.xplane;
			nds[9]=nodes[0]-mpmgrid.xplane;
			ndIDs[4]=4;
			ndIDs[5]=5;
			ndIDs[6]=6;
			ndIDs[7]=14;
			ndIDs[8]=15;
			*numnds=9;
		}
		else if(xipos->x<=q2)
		{	nds[5]=nodes[0]-mpmgrid.yplane;
			nds[6]=nds[5]+mpmgrid.xplane;
			ndIDs[4]=5;
			ndIDs[5]=6;
			*numnds=6;
		}
		else
		{	nds[5]=nodes[0]-mpmgrid.yplane;
			nds[6]=nds[5]+mpmgrid.xplane;
			nds[7]=nds[6]+mpmgrid.xplane;
			nds[8]=nodes[1]+mpmgrid.xplane;
			nds[9]=nodes[2]+mpmgrid.xplane;
			ndIDs[4]=5;
			ndIDs[5]=6;
			ndIDs[6]=7;
			ndIDs[7]=8;
			ndIDs[8]=9;
			*numnds=9;
		}
	}
	
	// middle two y quadrants
	else if(xipos->y<=q2)
	{	if(xipos->x<q1)
		{	nds[5]=nodes[3]-mpmgrid.xplane;
			nds[6]=nodes[0]-mpmgrid.xplane;
			ndIDs[4]=14;
			ndIDs[5]=15;
			*numnds=6;
		}
		else if(xipos->x<=q2)
		{	*numnds=4;
		}
		else
		{	nds[5]=nodes[1]+mpmgrid.xplane;
			nds[6]=nodes[2]+mpmgrid.xplane;
			ndIDs[4]=8;
			ndIDs[5]=9;
			*numnds=6;
		}
	}
	
	// upper y quadrant
	else
	{	if(xipos->x<q1)
		{	nds[5]=nodes[2]+mpmgrid.yplane;
			nds[6]=nds[5]-mpmgrid.xplane;
			nds[7]=nds[6]-mpmgrid.xplane;
			nds[8]=nodes[3]-mpmgrid.xplane;
			nds[9]=nodes[0]-mpmgrid.xplane;
			ndIDs[4]=11;
			ndIDs[5]=12;
			ndIDs[6]=13;
			ndIDs[7]=14;
			ndIDs[8]=15;
			*numnds=9;
		}
		else if(xipos->x<=q2)
		{	nds[5]=nodes[2]+mpmgrid.yplane;
			nds[6]=nds[5]-mpmgrid.xplane;
			ndIDs[4]=11;
			ndIDs[5]=12;
			*numnds=6;
		}
		else
		{	nds[5]=nodes[1]+mpmgrid.xplane;
			nds[6]=nodes[2]+mpmgrid.xplane;
			nds[7]=nodes[2]+mpmgrid.yplane+mpmgrid.xplane;
			nds[8]=nds[7]-mpmgrid.xplane;
			nds[9]=nds[8]-mpmgrid.xplane;
			ndIDs[4]=8;
			ndIDs[5]=9;
			ndIDs[6]=10;
			ndIDs[7]=11;
			ndIDs[8]=12;
			*numnds=9;
		}
	}
}

// get GIMP shape functions and optionally derivatives wrt x and y
// assumed to be properly numbered regular array
// input *xi position in element coordinate and ndIDs[0]... is which nodes (0-15)
void FourNodeIsoparam::GimpShapeFunction(Vector *xi,int numnds,int *ndIDs,int getDeriv,double *sfxn,
					double *xDeriv,double *yDeriv,double *zDeriv)
{
	int i;
	double xp,yp,Svpx,Svpy,dSvpx,dSvpy,xsign,ysign,argx=0.,argy=0.;
	
	// L is the cell spacing, 2*lp is the current particle size.
	// assuming the particle size is the same in x and y direction in element coordinate
	// the deformation of the particle is not considered yet.
	// assumes 4 particles per element
	double q1=0.5,q2=1.5,q3=2.5;
	
	for(i=0;i<numnds;i++)
	{	xp=fabs(xi->x-gxii[ndIDs[i]]);			// first qudrant (xp, yp)>=0
		yp=fabs(xi->y-geti[ndIDs[i]]);
		
		if(xp<=q1)
			Svpx=(4.*xp*xp-7.)/8.;
		else if(xp<=q2)
			Svpx=(xp-2.)/2.;
		else if(xp<=q3)
		{	argx=(5.-2.*xp)/4.;
			Svpx=-argx*argx;
		}
		else
			Svpx=0.;
			
		if(yp<=q1)
			Svpy=(4.*yp*yp-7.)/8.;
		else if(yp<=q2)
			Svpy=(yp-2.)/2.;
		else if(yp<=q3)
		{	argy=(5.-2.*yp)/4.;
			Svpy=-argy*argy;
		}
		else
			Svpy=0.;
			
		sfxn[i]=Svpx*Svpy;
				
		// find shape function at (xp,yp) 		
		if(getDeriv)
		{	xsign = xi->x>gxii[ndIDs[i]] ? 1. : -1.;
			ysign = xi->y>geti[ndIDs[i]] ? 1. : -1.;

			if(xp<=q1)
				dSvpx=xp;
			else if(xp<=q2)
				dSvpx=0.5;
			else if(xp<=q3)
				dSvpx=argx;
			else
				dSvpx=0.;
				
			if(yp<=q1)
				dSvpy=yp;
			else if(yp<=q2)
				dSvpy=0.5;
			else if(yp<=q3)
				dSvpy=argy;
			else
				dSvpy=0.;

			xDeriv[i]=xsign*dSvpx*Svpy*2.0/GetDeltaX();
			yDeriv[i]=ysign*Svpx*dSvpy*2.0/GetDeltaY();
		}
	}
}

#endif

#pragma mark FourNodeIsoparam::Accessors

// element name as an ID
short FourNodeIsoparam::ElementName(void) { return(FOUR_NODE_ISO); }

// number of nodes in this element
int FourNodeIsoparam::NumberNodes(void) { return 4; }
