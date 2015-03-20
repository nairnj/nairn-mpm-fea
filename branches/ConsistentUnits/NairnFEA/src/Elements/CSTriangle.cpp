/******************************************************************************** 
    CSTriangle.cpp
    NairnFEA
    
    Created by John Nairn on 10/24/05.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Elements/CSTriangle.hpp"
#include "NairnFEA_Class/NairnFEA.hpp"
#include "Materials/MaterialBase.hpp"

/********************************************************************************
	CSTriangle: Constructors and Destructor
********************************************************************************/

/* Main FEA constructor - passes on to Linear 2D
	But also sets extra node to match starting node
*/
CSTriangle::CSTriangle(int eNum,int *eNode,int eMat,double eAng,double eThick) : 
            Linear2D(eNum,eNode,eMat,eAng,eThick)
{
    nodes[3]=eNode[0];			// needed in some algorithms
}

/********************************************************************************
	CSTriangle: methods
********************************************************************************/

// element name as an ID
short CSTriangle::ElementName(void) { return CS_TRIANGLE; }

// number of nodes in this element
int CSTriangle::NumberNodes(void) const { return 3; }

// number of sides in this element
int CSTriangle::NumberSides(void) const { return 3; }

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
void CSTriangle::ShapeFunction(Vector *xi,int getDeriv,
		double *sfxn,double *xiDeriv,double *etaDeriv,Vector *eNodes,
                double *outDetjac,double *outAsr,double *asbe) const
{
	double ai,aj,ak,bi,bj,bk,ci,cj,ck,twicearea;
	
	ai=eNodes[1].x*eNodes[2].y-eNodes[2].x*eNodes[1].y;
	aj=eNodes[2].x*eNodes[0].y-eNodes[0].x*eNodes[2].y;
	ak=eNodes[0].x*eNodes[1].y-eNodes[1].x*eNodes[0].y;
	twicearea=(ai+aj+ak);
	bi=eNodes[1].y-eNodes[2].y;
	bj=eNodes[2].y-eNodes[0].y;
	bk=eNodes[0].y-eNodes[1].y;
	ci=eNodes[2].x-eNodes[1].x;
	cj=eNodes[0].x-eNodes[2].x;
	ck=eNodes[1].x-eNodes[0].x;
	sfxn[0]=(ai+bi*xi->x+ci*xi->y)/twicearea;
	sfxn[1]=(aj+bj*xi->x+cj*xi->y)/twicearea;
	sfxn[2]=(ak+bk*xi->x+ck*xi->y)/twicearea;
	if(getDeriv==BMATRIX)
	{	// elements of B Matrix
		xiDeriv[0]=bi/twicearea;
		xiDeriv[1]=bj/twicearea;
		xiDeriv[2]=bk/twicearea;
		etaDeriv[0]=ci/twicearea;
		etaDeriv[1]=cj/twicearea;
		etaDeriv[2]=ck/twicearea;
		
		// Ni/r if r≠0
		if(xi->x!=0.)
		{	asbe[0]=sfxn[0]/xi->x;
			asbe[1]=sfxn[1]/xi->x;
			asbe[2]=sfxn[2]/xi->x;
		}
		else
		{	asbe[0]=0.;
			asbe[1]=0.;
			asbe[2]=0.;
		}
		
		if(outDetjac!=NULL) *outDetjac=twicearea/2.;		/* Area of the triangle */
		if(outAsr!=NULL) *outAsr=(eNodes[0].x+eNodes[1].x+eNodes[2].x)/3.;	/* avg r */
	}
}

/* Calculate Stiffness Matrix
*/
void CSTriangle::Stiffness(int np)
{
	double detjac,asr,dv;
	double xiDeriv[MaxElNd],etaDeriv[MaxElNd],asbe[MaxElNd],sfxn[MaxElNd];
	double bte[MxFree*MaxElNd][5],temp;
	double thck=thickness,deltaT;
	int numnds=NumberNodes();
	int ind1,ind2,i,j,irow,jcol,nst=2*numnds;
	MaterialBase *matl=theMaterials[material-1];
	Vector xi;
	
    // Load nodal coordinates (ce[]), temperature (te[]), and
	//    material props (pr.C[][] and pr.alpha[])
    GetProperties(np);
    
    // Zero upper hald element stiffness matrix (se[]) and reaction vector (re[])
    for(irow=1;irow<=nst;irow++)
    {   re[irow]=0.;
        for(jcol=irow;jcol<=nst;jcol++)
            se[irow][jcol]=0.;
    }

	/* Call shape routine to calculate element B (be,asbe) matrix
		and the determinant of the Jacobian - both at centriod */
	xi.x=(ce[1].x+ce[2].x+ce[3].x)/3.;
	xi.y=(ce[1].y+ce[2].y+ce[3].y)/3.;
	ShapeFunction(&xi,BMATRIX,&sfxn[1],&xiDeriv[1],&etaDeriv[1],&ce[1],
		&detjac,&asr,&asbe[1]);
	
	/* Form matrix product BT E - exploit known sparcity of
		B matrix and only include multiplications by nonzero elements */
	deltaT=0.;
	ind1=-1;
	if(np!=AXI_SYM)
	{	dv=thck*detjac;
		for(i=1;i<=numnds;i++)
		{	ind1=ind1+2;
			ind2=ind1+1;
			for(j=1;j<=3;j++)
			{	bte[ind1][j]=dv*(xiDeriv[i]*matl->pr.C[1][j]+etaDeriv[i]*matl->pr.C[3][j]);
				bte[ind2][j]=dv*(etaDeriv[i]*matl->pr.C[2][j]+xiDeriv[i]*matl->pr.C[3][j]);
			}
			deltaT+=te[i]*sfxn[i];
		}
	}
	else
	{	dv=asr*detjac;
		for(i=1;i<=numnds;i++)
		{	ind1=ind1+2;
			ind2=ind1+1;
			for(j=1;j<=4;j++)
			{	bte[ind1][j]=dv*(xiDeriv[i]*matl->pr.C[1][j]+etaDeriv[i]*matl->pr.C[3][j]
								+asbe[i]*matl->pr.C[4][j]);
				bte[ind2][j]=dv*(etaDeriv[i]*matl->pr.C[2][j]+xiDeriv[i]*matl->pr.C[3][j]);
			}
			deltaT+=te[i]*sfxn[i];
		}
	}

	/* Form stiffness matrix by getting BT E B and add in initial
		strains into element load vector */
	for(irow=1;irow<=nst;irow++)
	{	for(j=1;j<=3;j++)
		{	re[irow]+=bte[irow][j]*matl->pr.alpha[j]*deltaT;
		}
        
		if(np!=AXI_SYM)
		{	for(jcol=irow;jcol<=nst;jcol++)
			{	if(IsEven(jcol))
				{	ind1=jcol/2;
					temp=bte[irow][2]*etaDeriv[ind1]
							+bte[irow][3]*xiDeriv[ind1];
				}
				else
				{	ind1=(jcol+1)/2;
					temp=bte[irow][1]*xiDeriv[ind1]
							+bte[irow][3]*etaDeriv[ind1];
				}
				se[irow][jcol]+=temp;
			}
		}
		else
		{	re[irow]+=bte[irow][4]*matl->pr.alpha[j]*deltaT;
			for(jcol=irow;jcol<=nst;jcol++)
			{	if(IsEven(jcol))
				{	ind1=jcol/2;
					temp=bte[irow][2]*etaDeriv[ind1]
							+bte[irow][3]*xiDeriv[ind1];
				}
				else
				{	ind1=(jcol+1)/2;
					temp=+bte[irow][1]*xiDeriv[ind1]
							+bte[irow][3]*etaDeriv[ind1]
							+bte[irow][4]*asbe[ind1];
				}
				se[irow][jcol]+=temp;
			}
		}
	}

	/* Fill in lower half of stiffness matrix */
	for(irow=1;irow<=nst-1;irow++)
	{	for(jcol=irow+1;jcol<=nst;jcol++)
		{	se[jcol][irow]=se[irow][jcol];
		}
	}
}

/* Calculate Element forces, stresses, and strain energy
*/
void CSTriangle::ForceStress(double *rm,int np,int nfree)
{
	double sgp[5],etot[5];
	double temp,dv;
	double xiDeriv[MaxElNd],etaDeriv[MaxElNd],asbe[MaxElNd],sfxn[MaxElNd];
	double detjac,asr;
	double thck=thickness,deltaT;
	int numnds=3,nst=2*numnds;
	int i,j,ind1,ind2,ind,indg;
    MaterialBase *matl=theMaterials[material-1];
	Vector xi;
    
    // Load element coordinates (ce[]), noodal temperature (te[]), 
	//    and material props (pr.C[][] and pr.alpha[])
    GetProperties(np);
    
    // Load nodal displacements into re[]
    ind=0;
    for(j=1;j<=numnds;j++)
    {	indg=nfree*(nodes[j-1]-1);
        for(i=1;i<=nfree;i++)
            re[++ind]=rm[indg+i];
    }

    // zero force at each degree of freedom (stored in se[i][7])
    for(i=1;i<=nst;i++)
        se[i][7]=0.;
        
    // element strainEnergy
    strainEnergy=0.;

	// Zero stress and force vectors to hold results
	nst=2*numnds;
	for(i=1;i<=4;i++) sgp[i]=0.;

	/* Call shape routine to calculate element B (be,asbe) matrix
		and the determinant of the Jacobian - both at centriod */
	xi.x=(ce[1].x+ce[2].x+ce[3].x)/3.;
	xi.y=(ce[1].y+ce[2].y+ce[3].y)/3.;
	ShapeFunction(&xi,BMATRIX,&sfxn[1],&xiDeriv[1],&etaDeriv[1],&ce[1],
		&detjac,&asr,&asbe[1]);
			
	// Evaluate volume element
	if(np!=AXI_SYM)
        dv=thck*detjac;
	else
        dv=asr*detjac;

	/* Evaluate etot=(B d - e0). In forming
		B d, avoid multiplications by zero */
	deltaT=0.;
	for(i=1;i<=numnds;i++) deltaT+=te[i]*sfxn[i];
	etot[1]=-matl->pr.alpha[1]*deltaT;
	etot[2]=-matl->pr.alpha[2]*deltaT;
	etot[3]=-matl->pr.alpha[3]*deltaT;
	etot[4]=-matl->pr.alpha[4]*deltaT;
	ind1=-1;
	for(i=1;i<=numnds;i++)
	{	ind1=ind1+2;
		ind2=ind1+1;
		etot[1]+=xiDeriv[i]*re[ind1];
		etot[2]+=etaDeriv[i]*re[ind2];
		etot[3]+=etaDeriv[i]*re[ind1]+xiDeriv[i]*re[ind2];
		if(np==AXI_SYM)
			etot[4]+=asbe[i]*re[ind1];
	}

	// Multply by stiffness matrix: sig = mdm * etot
	if(np!=AXI_SYM)
	{	for(i=1;i<=3;i++)
		{	for(j=1;j<=3;j++)
				sgp[i]+=matl->pr.C[i][j]*etot[j];
		}
	}
	else
	{	for(i=1;i<=4;i++)
		{	for(j=1;j<=4;j++)
				sgp[i]+=matl->pr.C[i][j]*etot[j];
		}
	}

	// Add terms for getting integral(BT sig) which gives nodal forces
	ind1=-1;
	for(i=1;i<=numnds;i++)
	{	ind1=ind1+2;
		ind2=ind1+1;
		temp=xiDeriv[i]*sgp[1]+etaDeriv[i]*sgp[3];
		if(np==AXI_SYM)
			temp+=asbe[i]*sgp[4];
		se[ind1][7]+=temp*dv;
		temp=etaDeriv[i]*sgp[2]+xiDeriv[i]*sgp[3];
		se[ind2][7]+=temp*dv;
	}

	// Get initial/thermal strain contribution to strain energy
	temp=0.;
	for(i=1;i<=3;i++) temp+=sgp[i]*matl->pr.alpha[i]*deltaT;
	if(np==AXI_SYM) temp+=sgp[4]*matl->pr.alpha[4]*deltaT;
	strainEnergy-=0.5*temp*dv;

	/* When plane strain account for constrained 1D shrinkage effect
			on strain energy */
	if(np==PLANE_STRAIN)
		strainEnergy+=0.5*matl->pr.C[4][4]*deltaT*deltaT*dv;
	
	// Add 1/2 Fd to strain energy
	temp=0.;
	for(i=1;i<=nst;i++) temp+=re[i]*se[i][7];
	strainEnergy+=0.5*temp;

	// copy stress to all nodes
	for(i=1;i<=numnds;i++)
	{	for(j=1;j<=4;j++)
			se[i][j]=sgp[j];
	}

	/* For plane strain analysis, calculate sigz stress
		For axisymmetric, multiply force and energy by 2 pi */
	if(np==PLANE_STRAIN)
    {   for(i=1;i<=numnds;i++)
			se[i][4]=matl->GetStressStrainZZ(se[i][1],se[i][2],se[i][3],te[i],angle,np);
    }
	else if(np==AXI_SYM)
	{	for(i=1;i<=nst;i++)
			se[i][7]*=2.*PI_CONSTANT;
		strainEnergy*=2*PI_CONSTANT;
	}
}

