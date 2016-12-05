/******************************************************************************** 
    Interface2D.cpp
    NairnFEA
    
    Created by John Nairn on 01/07/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Elements/Interface2D.hpp"
#include "NairnFEA_Class/NairnFEA.hpp"
#include "Materials/MaterialBase.hpp"

#pragma mark Interface2D: Constructors and Destructor

/* Main FEA constructor - passes on to Linear 2D
	But also sets extra node to match starting node
*/
Interface2D::Interface2D(int eNum,int *eNode,int eMat,double eAng,double eThick) : 
            ElementBase(eNum,eNode,eMat,eAng)
{
	thickness=eThick;
}

#pragma mark Interface2D: methods

//	Point can not be in an interface elemenbt
short Interface2D::PtInElement(Vector &pt) const { return FALSE; }

/* Calculate Element forces, stresses, and strain energy
*/
void Interface2D::ForceStress(double *rm,int np,int nfree)
{
	int numnds=NumberNodes();
	int i,j,nst=2*numnds;
    int indg;
	int nameEl=ElementName();
	double Fxy[(MaxElNd-1)*MxFree+1];
	double dx,dy,len,Dn,Dt;
	double pi=3.141592653589793,dsx,dsy,xpxi,ypxi;
    MaterialBase *matl=theMaterials[material-1];
	
	// Get stiffness matrix (and also load nodal coordinates (ce[]) and material props (pr.C[][]))
	Stiffness(np);
	
    // Load nodal displacements into re[]
    int ind=0;
    for(j=1;j<=numnds;j++)
    {	indg=nfree*(nodes[j-1]-1);
        for(i=1;i<=nfree;i++)
            re[++ind]=rm[indg+i];
    }

	// forces in cartensian frame
	for(i=1;i<=nst;i++)
	{	Fxy[i]=0.;
		for(j=1;j<=nst;j++)
		{	if(np!=AXI_SYM)
				Fxy[i]+=se[i][j]*re[j];
			else
				Fxy[i]+=2.*pi*se[i][j]*re[j];
		}
	}
	
	// transfer force to se[][] and get strain energy of interphase
    strainEnergy=0.;
	for(i=1;i<=nst;i++)
	{	strainEnergy+=0.5*re[i]*Fxy[i];
		se[i][7]=Fxy[i];
	}
	for(i=1;i<=numnds;i++)
	{	se[i][2]=0.;
		se[i][4]=0.;
	}
		
	// get normal and shear stresses from interface law
	Dn=matl->pr.C[1][1];
	Dt=matl->pr.C[1][2];
	if(nameEl==LINEAR_INTERFACE)
	{	dx=ce[2].x-ce[1].x;
		dy=ce[2].y-ce[1].y;
		len=sqrt(dx*dx + dy*dy);
		
		// nodes 1 and 4 and then 2 and 3
		InterfaceTraction(1,4,dx,dy,len,Dn,Dt);
		InterfaceTraction(2,3,dx,dy,len,Dn,Dt);
	}
	
	else
	{	dx=(ce[3].x-ce[1].x)/2.;
		dy=(ce[3].y-ce[1].y)/2.;
		dsx=ce[1].x+ce[3].x-2.*ce[2].x;
		dsy=ce[1].y+ce[3].y-2.*ce[2].y;
		
		// nodes 1 and 6
		xpxi=dx-dsx;
		ypxi=dy-dsy;
		len=sqrt(xpxi*xpxi + ypxi*ypxi);
		InterfaceTraction(1,6,xpxi,ypxi,len,Dn,Dt);
		
		// nodes 2 and 5
		xpxi=dx;
		ypxi=dy;
		len=sqrt(xpxi*xpxi + ypxi*ypxi);
		InterfaceTraction(2,5,xpxi,ypxi,len,Dn,Dt);
		
		// nodes 3 and 4
		xpxi=dx+dsx;
		ypxi=dy+dsy;
		len=sqrt(xpxi*xpxi + ypxi*ypxi);
		InterfaceTraction(3,4,xpxi,ypxi,len,Dn,Dt);
	}
}

// Calculate tractions at a node in interface element
void Interface2D::InterfaceTraction(int node1,int node2,double xpxi,double ypxi,double len,double Dn,double Dt)
{
	int y1=2*node1,x1=y1-1;
	int y2=2*node2,x2=y2-1;
	
	double un=(re[x1]-re[x2])*ypxi/len - (re[y1]-re[y2])*xpxi/len;
	double ut=(re[x1]-re[x2])*xpxi/len + (re[y1]-re[y2])*ypxi/len;
	double sig=Dn*un;
	double tau=Dt*ut;
	se[node1][1]=sig;
	se[node2][1]=sig;
	se[node1][3]=tau;
	se[node2][3]=tau;
}

// increment element of an interface stiffness matrix
void Interface2D::IncrementStiffnessElements(double dStiff,double *fn,
				double xpxi,double ypxi,double dlxi,double Dn,double Dt)
{
	int irow,jcol,nst=2*NumberNodes();
	int n1,n2,term;
	
	// loop over upper half diagonal of stiffness matrix
	for(irow=1;irow<=nst;irow++)
	{	for(jcol=irow;jcol<=nst;jcol++)
		{	// get term type and shape function numbers
			if(IsEven(irow))
			{	n1=irow/2;
				if(IsEven(jcol))
				{	n2=jcol/2;
					term=3;		// both even
				}
				else
				{	n2=(jcol+1)/2;
					term=2;		// even, odd
				}
			}
			else
			{	n1=(irow+1)/2;
				if(IsEven(jcol))
				{	n2=jcol/2;
					term=2;		// odd, even
				}
				else
				{	n2=(jcol+1)/2;
					term=1;		// both odd
				}
			}
			
			/* increment element depending on term  and current weight type
				Note: signs attached to shape functions take care of minus signs
						for matrix elements with displacement on opposite surface  */
			switch(term)
			{	case 1:			// both odd
					se[irow][jcol]+=dStiff*fn[n1]*fn[n2]*(Dt*xpxi*xpxi + Dn*ypxi*ypxi)/dlxi;
					break;
				case 2:			// one even and one odd
					se[irow][jcol]+=dStiff*fn[n1]*fn[n2]*(Dt-Dn)*xpxi*ypxi/dlxi;
					break;
				case 3:			// both even
					se[irow][jcol]+=dStiff*fn[n1]*fn[n2]*(Dt*ypxi*ypxi + Dn*xpxi*xpxi)/dlxi;
					break;
				default:
					break;
			}
		}
	}
}

//	Find extent of this element - but do not contribute to gridTolerance
void Interface2D::FindExtent(void)
{
    double saveGridTolerance = gridTolerance;
	ElementBase::FindExtent();
    gridTolerance = saveGridTolerance;
}

#pragma mark Interface2D: accessors

// number of sides in this element
int Interface2D::NumberSides(void) const { return 2; }

// thickness which may be in a subclass
double Interface2D::GetThickness(void) const { return thickness; }
void Interface2D::SetThickness(double thick) { thickness = thick; }

//	Interface element area is zero
double Interface2D::GetArea(void) const { return 0.; }

//	Interface element area is zero
double Interface2D::GetVolume(void) const { return 0.; }

// Bulk FEA element (some element may override with no)
bool Interface2D::BulkElement(void) { return FALSE; }

