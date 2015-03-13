/********************************************************************************
    EightNodeIsoparam.cpp
    NairnFEA
    
    Created by John Nairn on Fri Oct 22 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Elements/EightNodeIsoparam.hpp"

// Local globals
static double xii[4]={-1.,1.,1.,-1.};
static double eti[4]={-1.,-1.,1.,1.};

#pragma mark EightNodeIsoparam: Constructors and Destructor

/* Main FEA constructor - passes on to Quad 2D
	But also sets extra node to match starting node
*/
EightNodeIsoparam::EightNodeIsoparam(int eNum,int *eNode,int eMat,double eAng,double eThick) : 
            Quad2D(eNum,eNode,eMat,eAng,eThick)
{
    nodes[8]=eNode[0];			// needed in some algorithms
}

#pragma mark EightNodeIsoparam: methods

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
void EightNodeIsoparam::ShapeFunction(Vector *xi,int getDeriv,
		double *sfxn,double *xiDeriv,double *etaDeriv,Vector *eNodes,
                double *outDetjac,double *outAsr,double *asbe) const
{
	double temp1,temp2,asr;
	double jac[3][3],detjac;
	register int i;
	int ind1,ind2;
	
	for(i=0;i<4;i++)
	{	temp1=(1.+xii[i]*xi->x)/4.;
		temp2=(1.+eti[i]*xi->y)/4.;
		sfxn[i]=4.*temp1*temp2;
		if(getDeriv)
		{	xiDeriv[i]=xii[i]*temp2;
			etaDeriv[i]=eti[i]*temp1;
		}
	}
	sfxn[4]=(1.-xi->x*xi->x)*(1.-xi->y)/2.;
	sfxn[5]=(1.+xi->x)*(1.-xi->y*xi->y)/2.;
	sfxn[6]=(1.-xi->x*xi->x)*(1.+xi->y)/2.;
	sfxn[7]=(1.-xi->x)*(1.-xi->y*xi->y)/2.;
	if(getDeriv)
	{	xiDeriv[4]=-xi->x*(1.-xi->y);
		etaDeriv[4]=-(1.-xi->x*xi->x)/2.;
		xiDeriv[5]=(1.-xi->y*xi->y)/2.;
		etaDeriv[5]=-xi->y*(1.+xi->x);
		xiDeriv[6]=-xi->x*(1.+xi->y);
		etaDeriv[6]=(1.-xi->x*xi->x)/2.;
		xiDeriv[7]=-(1.-xi->y*xi->y)/2.;
		etaDeriv[7]=-xi->y*(1.-xi->x);
	}
	for(i=0;i<4;i++)
	{	if(i==0)
			ind1=7;
		else
			ind1=i+3;
		ind2=i+4;
		sfxn[i]-=(sfxn[ind1]+sfxn[ind2])/2.;
		if(getDeriv)
		{	xiDeriv[i]-=(xiDeriv[ind1]+xiDeriv[ind2])/2.;
			etaDeriv[i]-=(etaDeriv[ind1]+etaDeriv[ind2])/2.;
		}
	}
	
	// Get B matrix elements
	if(getDeriv==BMATRIX)
	{	// Find Jacobian Matrix
		jac[1][1]=0.;
		jac[1][2]=0.;
		jac[2][1]=0.;
		jac[2][2]=0.;
		for(i=0;i<8;i++)
		{	jac[1][1]+=xiDeriv[i]*eNodes[i].x;
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
		
		/* For radial position for axisymmetric analyses */
		asr=0.;
		for(i=0;i<8;i++)
			asr+=sfxn[i]*eNodes[i].x;
			
		/* Load J(-1)(1,1) ∂Ni/∂xi + J(-1)(1,2) ∂Ni/∂eta into xiDeriv[i]
		        J(-1)(2,1) ∂Ni/∂xi + J(-1)(2,2) ∂Ni/∂eta into etaDeriv[i]
		        Ni/r into asbe[i] */
		for(i=0;i<8;i++)
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

// Take stress at four gauss points and map them to 8 nodes in this element
// by using coordinate system on the gauss points (numbered ccw as 1,3,4,2)
// is mapping ot -1 to 1 coordinates.
// See FEA notes on stress extrapolation
// sgp[i][j] is stress j (1 to 4) at Gauss point i (1 to numGauss)
// se[i][j] is output stress j (1 to 4) at node i (1 to numnds) (externed variable)
void EightNodeIsoparam::ExtrapolateGaussStressToNodes(double sgp[][5])
{
	double gpt = 0.577350269189626;
	double at=1.+1./gpt;
	double bt=1.-1./gpt;
	double at2=at*at/4.;
	double bt2=bt*bt/4.;
	double atbt=at*bt/4.;
	
	double qe[9][5];
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
	
	at=at/4.;
	bt=bt/4.;
	qe[5][1]=at;
	qe[5][2]=bt;
	qe[5][3]=at;
	qe[5][4]=bt;
	qe[6][1]=bt;
	qe[6][2]=bt;
	qe[6][3]=at;
	qe[6][4]=at;
	qe[7][1]=bt;
	qe[7][2]=at;
	qe[7][3]=bt;
	qe[7][4]=at;
	qe[8][1]=at;
	qe[8][2]=at;
	qe[8][3]=bt;
	qe[8][4]=bt;
	
	// hard coded to 8 nodes and 4 Gauss points
	double temp;
	int i,j,k;
	for(i=1;i<=8;i++)
	{	for(j=1;j<=4;j++)
		{   temp=0.;
			for(k=1;k<=4;k++)
				temp+=qe[i][k]*sgp[k][j];
			se[i][j]=temp;
		}
	}
}

#pragma mark EightNodeIsoparam: accessors

// element name as an ID
short EightNodeIsoparam::ElementName(void) { return(EIGHT_NODE_ISO); }

// number of nodes in this element
int EightNodeIsoparam::NumberNodes(void) const { return 8; }

