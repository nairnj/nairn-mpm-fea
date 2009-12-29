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
EightNodeIsoparam::EightNodeIsoparam(long eNum,long *eNode,int eMat,double eAng,double eThick) : 
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
                double *outDetjac,double *outAsr,double *asbe)
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

#pragma mark EightNodeIsoparam: accessors

// element name as an ID
short EightNodeIsoparam::ElementName(void) { return(EIGHT_NODE_ISO); }

// number of nodes in this element
int EightNodeIsoparam::NumberNodes(void) { return 8; }

