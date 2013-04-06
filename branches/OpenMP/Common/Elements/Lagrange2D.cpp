/********************************************************************************
	Lagrange2D.cpp
	NairnFEA

	Created by John Nairn on 7 Feb 2011.
	Copyright (c) 2011 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Elements/Lagrange2D.hpp"
#include "Nodes/NodalPoint.hpp"
#ifdef MPM_CODE
	#include "NairnMPM_Class/MeshInfo.hpp"
#endif

// local variables
static int useForNumberNodes = 9;

#ifdef FEA_CODE
static int gaussOrder[9] = {1,7,9,3,4,8,6,2,5};
#endif

#pragma mark Lagrange2D: Constructors and Destructor

#ifdef MPM_CODE

// Main MPM constructor passes onto Quad2D constructor
Lagrange2D::Lagrange2D(int eNum,int *eNode) : Quad2D(eNum,eNode) {}

#else

// Main FEA constructor - passes on to Quad2D
Lagrange2D::Lagrange2D(int eNum,int *eNode,int eMat,double eAng,double eThick) : 
		Quad2D(eNum,eNode,eMat,eAng,eThick)
{
	numGauss=9;
	gaussSet=2;

}

#endif

#pragma mark Lagrange2D: methods

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
void Lagrange2D::ShapeFunction(Vector *xi,int getDeriv,
									  double *sfxn,double *xiDeriv,double *etaDeriv,Vector *eNodes,
									  double *outDetjac,double *outAsr,double *asbe) const
{
	sfxn[0]=(xi->x-1.)*xi->x * (xi->y-1.)*xi->y/4.;		// (-1,-1)
	sfxn[1]=(xi->x+1.)*xi->x * (xi->y-1.)*xi->y/4.;		// (1,-1)
	sfxn[2]=(xi->x+1.)*xi->x * (xi->y+1.)*xi->y/4.;		// (1,1)
	sfxn[3]=(xi->x-1.)*xi->x * (xi->y+1.)*xi->y/4.;		// (-1,1)
	sfxn[4]=(1.-xi->x*xi->x) * (xi->y-1.)*xi->y/2.;		// (0,-1)
	sfxn[5]=(xi->x+1.)*xi->x * (1.-xi->y*xi->y)/2.;		// (1,0)
	sfxn[6]=(1.-xi->x*xi->x) * (xi->y+1.)*xi->y/2.;		// (0,1)
	sfxn[7]=(xi->x-1.)*xi->x * (1.-xi->y*xi->y)/2.;		// (-1,0)
	sfxn[8]=(1.-xi->x*xi->x) * (1.-xi->y*xi->y);		// (0,0)
	
	if(getDeriv)
	{	xiDeriv[0]=(2.*xi->x-1.) * (xi->y-1.)*xi->y/4.;		// (-1,-1)
		xiDeriv[1]=(2.*xi->x+1.) * (xi->y-1.)*xi->y/4.;		// (1,-1)
		xiDeriv[2]=(2.*xi->x+1.) * (xi->y+1.)*xi->y/4.;		// (1,1)
		xiDeriv[3]=(2.*xi->x-1.) * (xi->y+1.)*xi->y/4.;		// (-1,1)
		xiDeriv[4]=(-2.*xi->x) * (xi->y-1.)*xi->y/2.;		// (0,-1)
		xiDeriv[5]=(2.*xi->x+1.) * (1.-xi->y*xi->y)/2.;		// (1,0)
		xiDeriv[6]=(-2.*xi->x) * (xi->y+1.)*xi->y/2.;		// (0,1)
		xiDeriv[7]=(2.*xi->x-1.) * (1.-xi->y*xi->y)/2.;		// (-1,0)
		xiDeriv[8]=(-2.*xi->x) * (1.-xi->y*xi->y);			// (0,0)
		
		etaDeriv[0]=(xi->x-1.)*xi->x * (2.*xi->y-1.)/4.;	// (-1,-1)
		etaDeriv[1]=(xi->x+1.)*xi->x * (2.*xi->y-1.)/4.;	// (1,-1)
		etaDeriv[2]=(xi->x+1.)*xi->x * (2.*xi->y+1.)/4.;	// (1,1)
		etaDeriv[3]=(xi->x-1.)*xi->x * (2.*xi->y+1.)/4.;	// (-1,1)
		etaDeriv[4]=(1.-xi->x*xi->x) * (2.*xi->y-1.)/2.;	// (0,-1)
		etaDeriv[5]=(xi->x+1.)*xi->x * (-2.*xi->y)/2.;		// (1,0)
		etaDeriv[6]=(1.-xi->x*xi->x) * (2.*xi->y+1.)/2.;	// (0,1)
		etaDeriv[7]=(xi->x-1.)*xi->x * (-2.*xi->y)/2.;		// (-1,0)
		etaDeriv[8]=(1.-xi->x*xi->x) * (-2.*xi->y);			// (0,0)
	}
	
	// Get B matrix elements
	if(getDeriv==BMATRIX)
	{	double jac[3][3];
		int i;
		
		// Find Jacobian Matrix
		jac[1][1]=0.;
		jac[1][2]=0.;
		jac[2][1]=0.;
		jac[2][2]=0.;
		for(i=0;i<9;i++)
		{	jac[1][1]+=xiDeriv[i]*eNodes[i].x;
			jac[1][2]+=xiDeriv[i]*eNodes[i].y;
			jac[2][1]+=etaDeriv[i]*eNodes[i].x;
			jac[2][2]+=etaDeriv[i]*eNodes[i].y;
		}
		double detjac=jac[1][1]*jac[2][2]-jac[1][2]*jac[2][1];
		double temp1=jac[1][1]/detjac,temp2;
		jac[1][1]=jac[2][2]/detjac;
		jac[2][2]=temp1;
		jac[1][2]=-jac[1][2]/detjac;
		jac[2][1]=-jac[2][1]/detjac;
		
		/* For radial position for axisymmetric analyses */
		double asr=0.;
		for(i=0;i<9;i++)
			asr+=sfxn[i]*eNodes[i].x;
		
		/* Load J(-1)(1,1) ∂Ni/∂xi + J(-1)(1,2) ∂Ni/∂eta into xiDeriv[i]
		 J(-1)(2,1) ∂Ni/∂xi + J(-1)(2,2) ∂Ni/∂eta into etaDeriv[i]
		 Ni/r into asbe[i] */
		for(i=0;i<9;i++)
		{	temp1=xiDeriv[i];
			temp2=etaDeriv[i];
			xiDeriv[i]=jac[1][1]*temp1 + jac[1][2]*temp2;
			etaDeriv[i]=jac[2][1]*temp1 + jac[2][2]*temp2;
			if(asr!=0.)
				asbe[i]=sfxn[i]/asr;
			else
				asbe[i]=0.;
		}
		
		// save other results
		if(outDetjac!=NULL) *outDetjac=detjac;
		if(outAsr!=NULL) *outAsr=asr;
	}
}

#ifdef FEA_CODE

// Take stress at 9 gauss points and map them to 9 nodes in this element
// by using coordinate system on the gauss points (numbered as element as
//		1,7,9,3,4,8,6,2,5) as mapping to -1 to 1 coordinates.
// sgp[i][j] is stress j (1 to 4) at Gauss point i (1 to numGauss)
// se[i][j] is output stress j (1 to 4) at node i (1 to numnds) (externed variable)
void Lagrange2D::ExtrapolateGaussStressToNodes(double sgp[][5])
{
	double gpt = 1./0.7745966692414834;
	double temp,sfxn[10];
	int i,j,k;
	Vector xi;
	
	// loop over 9 nodes
	for(i=1;i<=9;i++)
	{	if(i==1 || i==8 || i==4)
			xi.x = -gpt;
		else if(i==5 || i==9 || i==7)
			xi.x = 0.;
		else
			xi.x = gpt;
		if(i==1 || i==5 || i==2)
			xi.y = -gpt;
		else if(i==8 || i==9 || i==6)
			xi.y = 0.;
		else
			xi.y = gpt;
		
		ShapeFunction(&xi,FALSE,sfxn,NULL,NULL,NULL,NULL,NULL,NULL);
		for(j=1;j<=4;j++)
		{	temp=0.;
			for(k=0;k<9;k++)
				temp+=sfxn[k]*sgp[gaussOrder[k]][j];
			se[i][j] = temp;
		}
	}
}

#endif

#ifdef MPM_CODE

// get shape functions and optionally derivatives wrt x and y, but derivatives only work
// if it is a regular array. Shape functions are general
void Lagrange2D::ShapeFunction(Vector *xi,int getDeriv,double *sfxn,
									 double *xDeriv,double *yDeriv,double *zDeriv) const
{
	sfxn[0]=0.25 * (xi->x-1.)*xi->x * (xi->y-1.)*xi->y;		// (-1,-1)
	sfxn[1]=0.25 * (xi->x+1.)*xi->x * (xi->y-1.)*xi->y;		// (1,-1)
	sfxn[2]=0.25 * (xi->x+1.)*xi->x * (xi->y+1.)*xi->y;		// (1,1)
	sfxn[3]=0.25 * (xi->x-1.)*xi->x * (xi->y+1.)*xi->y;		// (-1,1)
	sfxn[4]=0.5 * (1.-xi->x*xi->x) * (xi->y-1.)*xi->y;		// (0,-1)
	sfxn[5]=0.5 * (xi->x+1.)*xi->x * (1.-xi->y*xi->y);		// (1,0)
	sfxn[6]=0.5 * (1.-xi->x*xi->x) * (xi->y+1.)*xi->y;		// (0,1)
	sfxn[7]=0.5 * (xi->x-1.)*xi->x * (1.-xi->y*xi->y);		// (-1,0)
	sfxn[8]=(1.-xi->x*xi->x) * (1.-xi->y*xi->y);			// (0,0)
	
	if(getDeriv)
	{	xDeriv[0]=0.5 * (2.*xi->x-1.) * (xi->y-1.)*xi->y/GetDeltaX();		// (-1,-1)
		xDeriv[1]=0.5 * (2.*xi->x+1.) * (xi->y-1.)*xi->y/GetDeltaX();		// (1,-1)
		xDeriv[2]=0.5 * (2.*xi->x+1.) * (xi->y+1.)*xi->y/GetDeltaX();		// (1,1)
		xDeriv[3]=0.5 * (2.*xi->x-1.) * (xi->y+1.)*xi->y/GetDeltaX();		// (-1,1)
		xDeriv[4]=(-2.*xi->x) * (xi->y-1.)*xi->y/GetDeltaX();				// (0,-1)
		xDeriv[5]=(2.*xi->x+1.) * (1.-xi->y*xi->y)/GetDeltaX();			// (1,0)
		xDeriv[6]=(-2.*xi->x) * (xi->y+1.)*xi->y/GetDeltaX();				// (0,1)
		xDeriv[7]=(2.*xi->x-1.) * (1.-xi->y*xi->y)/GetDeltaX();			// (-1,0)
		xDeriv[8]=2.0 * (-2.*xi->x) * (1.-xi->y*xi->y)/GetDeltaX();		// (0,0)
		
		yDeriv[0]=0.5 * (xi->x-1.)*xi->x * (2.*xi->y-1.)/GetDeltaY();		// (-1,-1)
		yDeriv[1]=0.5 * (xi->x+1.)*xi->x * (2.*xi->y-1.)/GetDeltaY();		// (1,-1)
		yDeriv[2]=0.5 * (xi->x+1.)*xi->x * (2.*xi->y+1.)/GetDeltaY();		// (1,1)
		yDeriv[3]=0.5 * (xi->x-1.)*xi->x * (2.*xi->y+1.)/GetDeltaY();		// (-1,1)
		yDeriv[4]=(1.-xi->x*xi->x) * (2.*xi->y-1.)/GetDeltaY();			// (0,-1)
		yDeriv[5]=(xi->x+1.)*xi->x * (-2.*xi->y)/GetDeltaY();				// (1,0)
		yDeriv[6]=(1.-xi->x*xi->x) * (2.*xi->y+1.)/GetDeltaY();			// (0,1)
		yDeriv[7]=(xi->x-1.)*xi->x * (-2.*xi->y)/GetDeltaY();				// (-1,0)
		yDeriv[8]=2.0 * (1.-xi->x*xi->x) * (-2.*xi->y)/GetDeltaY();		// (0,0)
	}
}

//	Find extent of this element - called once at start (and must be called)
void Lagrange2D::FindExtent(void)
{	
	ElementBase::FindExtent();
	
	/* for speed in GetXiPos() calculations - if it a parallelogram
	 Note: assumes element does not move. If it does, must recalculate these terms
	 */
	pgElement=TRUE;
	// are edges parallel wrt x coordinate?
	double xdel=fabs(nd[nodes[2]]->x-nd[nodes[1]]->x);
    double ydel=fabs(nd[nodes[3]]->x-nd[nodes[0]]->x);
	if(!DbleEqual(xdel,ydel) && xdel>1.e-10 && ydel>1.e-10)
	{	pgElement=FALSE;
		return;
	}

	// are edges parallel wrt y coordinate?
	xdel=fabs(nd[nodes[2]]->y-nd[nodes[1]]->y);
	ydel=fabs(nd[nodes[3]]->y-nd[nodes[0]]->y);
	if(!DbleEqual(xdel,ydel) && xdel>1.e-10 && ydel>1.e-10)
	{	pgElement=FALSE;
		return;
	}

	// mid side nodes in line?
	int i;
	for(i=0;i<3;i++)
	{	if(!DbleEqual(2.*nd[nodes[i+4]]->x, nd[nodes[i]]->x+nd[nodes[i+1]]->x) ||
				!DbleEqual(2.*nd[nodes[i+4]]->y, nd[nodes[i]]->y+nd[nodes[i+1]]->y))
		{	pgElement=FALSE;
			return;
		}
	}
	if(!DbleEqual(2.*nd[nodes[7]]->x, nd[nodes[3]]->x+nd[nodes[0]]->x) ||
		   !DbleEqual(2.*nd[nodes[7]]->y, nd[nodes[3]]->y+nd[nodes[0]]->y))
	{	pgElement=FALSE;
		return;
	}
		
	// Is center node in the center?
	if(!DbleEqual(2.*nd[nodes[8]]->x, nd[nodes[4]]->x+nd[nodes[6]]->x) ||
	   !DbleEqual(2.*nd[nodes[8]]->y, nd[nodes[4]]->y+nd[nodes[6]]->y))
	{	pgElement=FALSE;
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
void Lagrange2D::GetXiPos(Vector *pos,Vector *xipos) const
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
int Lagrange2D::Orthogonal(double *dx,double *dy,double *dz)
{
	int i;
	double xdel,ydel,xdelmid,ydelmid;
	
	*dx=0.;
	*dy=0.;
	*dz=0.;
	for(i=0;i<3;i++)
    {	xdelmid=fabs(nd[nodes[i+4]]->x-nd[nodes[i]]->x);
		ydelmid=fabs(nd[nodes[i+4]]->y-nd[nodes[i]]->y);
    	xdel=fabs(nd[nodes[i+1]]->x-nd[nodes[i]]->x);
		ydel=fabs(nd[nodes[i+1]]->y-nd[nodes[i]]->y);
		if((xdel>=1.e-12 && ydel>=1.e-12) || (xdelmid>=1.e-12 && ydelmid>=1.e-12)) return FALSE;
		*dx=fmax(xdel,*dx);
		*dy=fmax(ydel,*dy);
	}
	xdelmid=fabs(nd[nodes[7]]->x-nd[nodes[3]]->x);
	ydelmid=fabs(nd[nodes[7]]->y-nd[nodes[3]]->y);
    xdel=fabs(nd[nodes[0]]->x-nd[nodes[3]]->x);
	ydel=fabs(nd[nodes[0]]->y-nd[nodes[3]]->y);
	if((xdel>=1.e-12 && ydel>=1.e-12) || (xdelmid>=1.e-12 && ydelmid>=1.e-12)) return FALSE;
	*dx=fmax(xdel,*dx);
	*dy=fmax(ydel,*dy);
	return TRUE;
}

#endif

#pragma mark Lagrange2D: accessors

// trick into 8 node temporarily for area calculation
double Lagrange2D::GetArea(void) const
{	useForNumberNodes=8;
	double theArea = Quad2D::GetArea();
	useForNumberNodes=9;
	return theArea;
}

// element name as an ID
short Lagrange2D::ElementName(void) { return NINE_NODE_LAGRANGE; }

// number of nodes in this element
int Lagrange2D::NumberNodes(void) const { return useForNumberNodes; }
