/******************************************************************************** 
    SixNodeTriangle.cpp
    NairnFEA
    
    Created by John Nairn on Mon Oct 25 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Elements/SixNodeTriangle.hpp"

#pragma mark SixNodeTriangle: Constructors and Destructor

/* Main FEA constructor - passes on to Quad 2D
	But also sets extra node to match starting node
*/
SixNodeTriangle::SixNodeTriangle(int eNum,int *eNode,int eMat,double eAng,double eThick) : 
            Quad2D(eNum,eNode,eMat,eAng,eThick)
{
    nodes[6]=eNode[0];			// needed in some algorithms
	numGauss=3;
	gaussSet=1;
}

#pragma mark SixNodeTriangle: methods

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
void SixNodeTriangle::ShapeFunction(Vector *xi,int getDeriv,
		double *sfxn,double *xiDeriv,double *etaDeriv,Vector *eNodes,
                double *outDetjac,double *outAsr,double *asbe) const
{
	double jac[3][3],detjac,asr,temp1,temp2;
	int i;
	double pxi3;
	
	pxi3=1.-xi->x-xi->y;
	sfxn[0]=xi->x*(2.*xi->x-1.);
	sfxn[1]=xi->y*(2.*xi->y-1.);
	sfxn[2]=pxi3*(2.*pxi3-1.);
	sfxn[3]=4.*xi->x*xi->y;
	sfxn[4]=4.*xi->y*pxi3;
	sfxn[5]=4.*xi->x*pxi3;
	
	if(getDeriv || getDeriv==BMATRIX)
	{	xiDeriv[0]=4*xi->x-1.;
		xiDeriv[1]=0.;
		xiDeriv[2]=-4.*pxi3+1.;
		xiDeriv[3]=4.*xi->y;
		xiDeriv[4]=-4.*xi->y;
		xiDeriv[5]=4.*(pxi3-xi->x);
		etaDeriv[0]=0.;
		etaDeriv[1]=4*xi->y-1.;
		etaDeriv[2]=-4.*pxi3+1.;
		etaDeriv[3]=4.*xi->x;
		etaDeriv[4]=4.*(pxi3-xi->y);
		etaDeriv[5]=-4.*xi->x;
	}

	// Get B matrix elements
	if(getDeriv==BMATRIX)
	{	// Find Jacobian Matrix
		jac[1][1]=0.;
		jac[1][2]=0.;
		jac[2][1]=0.;
		jac[2][2]=0.;
		for(i=0;i<6;i++)
		{	jac[1][1]+=xiDeriv[i]*eNodes[i].x;
			jac[1][2]+=xiDeriv[i]*eNodes[i].y;
			jac[2][1]+=etaDeriv[i]*eNodes[i].x;
			jac[2][2]+=etaDeriv[i]*eNodes[i].y;
		}
		detjac=jac[1][1]*jac[2][2]-jac[1][2]*jac[2][1];

		// invert Jacobian
		temp1=jac[1][1]/detjac;
		jac[1][1]=jac[2][2]/detjac;
		jac[2][2]=temp1;
		jac[1][2]=-jac[1][2]/detjac;
		jac[2][1]=-jac[2][1]/detjac;
		
		// For radial position for axisymmetric analyses
		asr=0.;
		for(i=0;i<6;i++)
			asr=asr+sfxn[i]*eNodes[i].x;
			
		/* Load J(-1)(1,1) ∂Ni/∂xi + J(-1)(1,2) ∂Ni/∂eta into xiDeriv[i]
		        J(-1)(2,1) ∂Ni/∂xi + J(-1)(2,2) ∂Ni/∂eta into etaDeriv[i]
		        Ni/r into asbe[i] */
		for(i=0;i<6;i++)
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

// Take stress at three gauss points and map them to 6 nodes in this element
// See FEA notes on stress extrapolation
// sgp[i][j] is stress j (1 to 4) at Gauss point i (1 to numGauss)
// se[i][j] is output stress j (1 to 4) at node i (1 to numnds) (externed variable)
void SixNodeTriangle::ExtrapolateGaussStressToNodes(double sgp[][5])
{
	// extraplate internal triangle to 6 nodes - see notes FEA section
	int j;
	for(j=1;j<=4;j++)
	{	se[1][j]=(5.*sgp[1][j]-sgp[2][j]-sgp[3][j])/3.;
		se[2][j]=(-sgp[1][j]+5.*sgp[2][j]-sgp[3][j])/3.;
		se[3][j]=(-sgp[1][j]-sgp[2][j]+5.*sgp[3][j])/3.;
		se[4][j]=(2.*sgp[1][j]+2.*sgp[2][j]-sgp[3][j])/3.;
		se[5][j]=(-sgp[1][j]+2.*sgp[2][j]+2.*sgp[3][j])/3.;
		se[6][j]=(2.*sgp[1][j]-sgp[2][j]+2.*sgp[3][j])/3.;
	}
}

#pragma mark SixNodeTriangle: accessors

// element name as an ID
short SixNodeTriangle::ElementName(void) { return(ISO_TRIANGLE); }

// number of nodes in this element
int SixNodeTriangle::NumberNodes(void) const { return 6; }

// number of sides in this element (override if differs)
int SixNodeTriangle::NumberSides(void) const { return(3); }

