/********************************************************************************
    EdgeBC.cpp
    NairnFEA
    
    Created by John Nairn on Mon Mar 15 2004.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Boundary_Conditions/EdgeBC.hpp"
#include "Elements/ElementBase.hpp"

// Nodal BC global
EdgeBC *firstEdgeBC=NULL;

#pragma mark EdgeBC: Constructors and Destructors

EdgeBC::EdgeBC(int elemNum,int faceNum,int dirNum) : FEABoundaryCondition()
{
	element=elemNum;
	face=faceNum;
	direction=dirNum;
	stress[0]=stress[1]=stress[2]=0.;
}

#pragma mark EdgeBC: Methods

// print load BC
EdgeBC *EdgeBC::PrintEdgeLoad(void)
{
	char loadType[10];
    char nline[200];
    size_t nlsize=200;
	
	if(direction==1)
		strcpy(loadType,"Normal");
	else
		strcpy(loadType,"Shear1");
	
	snprintf(nline,nlsize,"%5d  %2d  %6s",element,face,loadType);
    cout << nline;
	
	switch(facenodes)
	{   case 2:
            snprintf(nline,nlsize,"  %15.7e  %15.7e",stress[0],stress[1]);
			break;
		case 3:
            snprintf(nline,nlsize,"  %15.7e  %15.7e  %15.7e",stress[0],stress[1],stress[2]);
			break;
		default:
			break;
	}
	cout << nline << endl;
	
	return (EdgeBC *)nextObject;
}

// find consistent loads for this boundary condition
void EdgeBC::GetConsistentLoads(double *re,int np)
{
	theElements[element-1]->CalcEdgeLoads(re,face,direction,stress,np);
}

#pragma mark EdgeBC: Accessors

// set stress - one for each face node
void EdgeBC::SetStress(double *sigma,int num)
{
	int i;
	for(i=0;i<num;i++) stress[i]=sigma[i];
	facenodes=num;
}

// equal means applied to same edge and dof, value does not matter
bool EdgeBC::SameDofSetting(EdgeBC *rhs)
{	return (element==rhs->element) && (face==rhs->face) && (direction==rhs->direction);
}

// zero based element index
int EdgeBC::ElementIndex(void) { return element-1; }

