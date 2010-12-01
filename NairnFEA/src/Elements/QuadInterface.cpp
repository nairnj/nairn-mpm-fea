/******************************************************************************** 
    QuadInterface.cpp
    NairnFEA
    
    Created by John Nairn on 01/07/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Elements/QuadInterface.hpp"
#include "NairnFEA_Class/NairnFEA.hpp"
#include "Materials/MaterialBase.hpp"

#define GAUSS_INT_PTS 4
static double gxi[4]={-0.861136311594053,-0.339981043584856,
				0.339981043584856,0.861136311594053};
static double gwt[4]={0.347854845137448,0.652145154861630,
				0.652145154861630, 0.347854845137448};

#pragma mark QuadInterface: Constructors and Destructor

/* Main FEA constructor - passes on to Linear 2D
	But also sets extra node to match starting node
*/
QuadInterface::QuadInterface(int eNum,int *eNode,int eMat,double eAng,double eThick) : 
            Interface2D(eNum,eNode,eMat,eAng,eThick)
{
    nodes[6]=eNode[0];			// needed in some algorithms
}

#pragma mark QuadInterface: methods

/* Calculate Shape Functions
*/
void QuadInterface::ShapeFunction(Vector *xi,int getDeriv,
		double *sfxn,double *xiDeriv,double *etaDeriv,Vector *eNodes,
                double *outDetjac,double *outAsr,double *asbe)
{
	// shape functions
	sfxn[0]=(xi->x*xi->x-xi->x)/2.;
	sfxn[1]=1.-xi->x*xi->x;
	sfxn[2]=(xi->x*xi->x+xi->x)/2.;
	sfxn[3]=-sfxn[2];
	sfxn[4]=-sfxn[1];
	sfxn[5]=-sfxn[0];
	
	// Find radial position for axisymmetric analyses
	if(outAsr!=NULL)
		*outAsr=sfxn[0]*eNodes[0].x + sfxn[1]*eNodes[1].x + sfxn[2]*eNodes[2].x;
}

/* Calculate Stiffness Matrix
*/
void QuadInterface::Stiffness(int np)
{
	double xpxi,ypxi,dlxi,fn[MaxElNd],temp,asr;
	Vector xi;
	
    // Load nodal coordinates (ce[]) and material props (mdm[][])
    GetProperties(np);
    
    // Zero stiffness (se[][]) and reaction (re[])
	ZeroUpperHalfStiffness();

	// basic parameters */
    MaterialBase *matl=theMaterials[material-1];
	double Dn=matl->mdm[1][1];
	double Dt=matl->mdm[1][2];
	double dx=(ce[3].x-ce[1].x)/2.;
	double dy=(ce[3].y-ce[1].y)/2.;
	double dsx=ce[1].x+ce[3].x-2.*ce[2].x;
	double dsy=ce[1].y+ce[3].y-2.*ce[2].y;
	
	int i;
	for(i=0;i<GAUSS_INT_PTS;i++)
	{	// get shape function at next xi
		xi.x=gxi[i];
		ShapeFunction(&xi,FALSE,&fn[1],NULL,NULL,&ce[1],NULL,&asr,NULL);
		
		// get integration terms that vary in xi
		xpxi=dx+dsx*xi.x;
		ypxi=dy+dsy*xi.x;
		dlxi=sqrt(xpxi*xpxi + ypxi*ypxi);
		
		// total weight
		if(np!=AXI_SYM)
			temp=gwt[i]*GetThickness()/1000.;
		else
			temp=gwt[i]*asr;		// stiffness matrix is force per radian, hence no 2π
		
		// increment all stiffness elements at this point
		IncrementStiffnessElements(temp,fn,xpxi,ypxi,dlxi,Dn,Dt);
	}
		
	// Fill in lower half of stiffness matrix
	FillLowerHalfStiffness();
}

#pragma mark QuadInterface: accessors

// element name as an ID
short QuadInterface::ElementName(void) { return QUAD_INTERFACE; }

// number of nodes in this element
int QuadInterface::NumberNodes(void) { return 6; }

// face nodes
int QuadInterface::FaceNodes(void) { return 3; }

