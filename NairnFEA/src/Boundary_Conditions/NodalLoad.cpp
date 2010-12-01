/********************************************************************************
    NodalLoad.cpp
    NairnFEA
    
    Created by John Nairn on Wed Mar 11 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Boundary_Conditions/NodalLoad.hpp"

// Nodal BC global
NodalLoad *firstLoadBC=NULL;

#pragma mark NodalLoad: Constructors and Destructors

NodalLoad::NodalLoad(int num,int dof) : FEABoundaryCondition()
{
    nodeNum=num;
    direction=dof;
}

#pragma mark NodalLoad: Methods

// print load BC
NodalLoad *NodalLoad::PrintLoad(void)
{
    char nline[200];
    
    sprintf(nline,"%5d  %2d  %15.7e",nodeNum,direction,bcValue);
    cout << nline << endl;
	
	return (NodalLoad *)nextObject;
}

// remap if resequenced
NodalLoad *NodalLoad::MapNodes(int *revMap)
{
	nodeNum=revMap[nodeNum];
	return (NodalLoad *)nextObject;
}

// put load into reaction vector
NodalLoad *NodalLoad::Reaction(double *rm,int np,int nfree)
{
    int nglob;
    
    nglob=nfree*(nodeNum-1)+direction;
    if(np!=AXI_SYM)
        rm[nglob]=bcValue;
    else
        rm[nglob]=bcValue/(2.*PI_CONSTANT);
	
	return (NodalLoad *)nextObject;
}

#pragma mark NodalLoad: Methods

// equal means applied to same node and dof, value does not matter
bool NodalLoad::SameDofSetting(NodalLoad *rhs)
{	return (nodeNum==rhs->nodeNum) && (direction==rhs->direction);
}

