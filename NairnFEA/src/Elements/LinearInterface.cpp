/******************************************************************************** 
    LinearInterface.cpp
    NairnFEA
    
    Created by John Nairn on 01/07/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Elements/LinearInterface.hpp"
#include "NairnFEA_Class/NairnFEA.hpp"
#include "Materials/MaterialBase.hpp"

#define GAUSS_INT_PTS 4
static double gxi[4]={-0.861136311594053,-0.339981043584856,
				0.339981043584856,0.861136311594053};
static double gwt[4]={0.347854845137448,0.652145154861630,
				0.652145154861630, 0.347854845137448};

#pragma mark LinearInterface: Constructors and Destructor

/* Main FEA constructor - passes on to Linear 2D
	But also sets extra node to match starting node
*/
LinearInterface::LinearInterface(long eNum,long *eNode,int eMat,double eAng,double eThick) : 
            Interface2D(eNum,eNode,eMat,eAng,eThick)
{
    nodes[4]=eNode[0];			// needed in some algorithms
}

#pragma mark LinearInterface: methods

/* Calculate Shape Functions
*/
void LinearInterface::ShapeFunction(Vector *xi,int getDeriv,
		double *sfxn,double *xiDeriv,double *etaDeriv,Vector *eNodes,
                double *outDetjac,double *outAsr,double *asbe)
{
	// shape functions
	sfxn[0]=(1.-xi->x)/2.;
	sfxn[1]=(1.+xi->x)/2.;
	sfxn[2]=-sfxn[1];
	sfxn[3]=-sfxn[0];
	
	// Finf radial position for axisymmetric analyses
	if(outAsr!=NULL)
		*outAsr=sfxn[0]*eNodes[0].x + sfxn[1]*eNodes[1].x;
}

/* Calculate Stiffness Matrix
*/
void Interface2D::Stiffness(int np)
{
	double fn[MaxElNd],temp,asr;
	Vector xi;
	
    // Load nodal coordinates (ce[]) and material props (mdm[][])
    GetProperties(np);
    
    // Zero stiffness (se[][]) and reaction (re[])
	ZeroUpperHalfStiffness();

	// basic parameters */
    MaterialBase *matl=theMaterials[material-1];
	double Dn=matl->mdm[1][1];
	double Dt=matl->mdm[1][2];
	double xpxi=(ce[2].x-ce[1].x)/2;
	double ypxi=(ce[2].y-ce[1].y)/2;
	double dlxi=sqrt(xpxi*xpxi + ypxi*ypxi);
	
	int i;
	for(i=0;i<GAUSS_INT_PTS;i++)
	{	// get shape function at next xi
		xi.x=gxi[i];
		ShapeFunction(&xi,FALSE,&fn[1],NULL,NULL,&ce[1],NULL,&asr,NULL);
		
		// total weight
		if(np!=AXI_SYM)
			temp=gwt[i]*GetThickness()/1000.;
		else
			temp=gwt[i]*asr;		// stiffness matrix is force per radian, hence no 2 pi
				
		// increment all stiffness elements at this point
		IncrementStiffnessElements(temp,fn,xpxi,ypxi,dlxi,Dn,Dt);
	}
		
	// Fill in lower half of stiffness matrix
	FillLowerHalfStiffness();
}

#pragma mark LinearInterface: accessors

// element name as an ID
short LinearInterface::ElementName(void) { return LINEAR_INTERFACE; }

// number of nodes in this element
int LinearInterface::NumberNodes(void) { return 4; }

// face nodes
int LinearInterface::FaceNodes(void) { return 2; }

