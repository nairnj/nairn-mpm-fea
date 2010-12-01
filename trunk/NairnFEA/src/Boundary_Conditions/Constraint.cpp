/********************************************************************************
    Constraint.cpp
    NairnFEA
    
    Created by John Nairn on 2/6/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Boundary_Conditions/Constraint.hpp"

// Constraint global
Constraint *firstConstraint=NULL;
int numConstraints=0;

#pragma mark Constraint: Constructors and Destructors

// baseNodes is contrained node, dofNum is 1 or 2 for x or y
// must set lambda number later
Constraint::Constraint(int baseNode,int dofNum)
{
	nodes.push_back(baseNode);
	coeff.push_back((double)1.);
	dof=dofNum;
	qvalue=0.;
}

#pragma mark Constraint: Methods

// Enter as right side of constrain equation, but add as -factor*node
void Constraint::AddNode(int node,double factor)
{
	nodes.push_back(node);
	coeff.push_back(-factor);
}

// set the constant
void Constraint::AddConstant(double factor) { qvalue=factor; }

// remap if resequenced
Constraint *Constraint::MapNodes(int *revMap)
{
	unsigned i;
	for(i=0;i<nodes.size();i++)
		nodes[i]=revMap[nodes[i]];
	return (Constraint *)nextObject;
}

// bandwidth calculation - nbase is number of dof not counting constraints
int Constraint::GetBandWidth(int nbase,int nfree)
{
	int cmax=nbase+lambda-nfree*(nodes[0]-1)-dof;
	
	unsigned i;
	for(i=1;i<nodes.size();i++)
		cmax=fmax(cmax,nbase+lambda-nfree*(nodes[i]-1)-dof);
	
	return cmax+1;		// one extra for diagonal
}

#pragma mark Constraint: Accessors

// number of noddes in the constraint
int Constraint::NumberNodes(void) { return (int)nodes.size(); }

// nodal DOF and coefficient (i is 1-based but stored here 0 based)
int Constraint::NodalDof(int i,int nfree) { return nfree*(nodes[i-1]-1)+dof; }
double Constraint::GetCoeff(int i) { return coeff[i-1]; }
double Constraint::GetQ(void) { return qvalue; }

// constraint number
void Constraint::SetLambdaNum(int num) { lambda=num; }
int Constraint::GetLambdaNum(void) { return lambda; }

// for debugging
void Constraint::Describe(void)
{	cout << "# Constraint " << lambda << " on u(" << dof << ") using nd[" << nodes[0] << "] = ";
	unsigned i;
	for(i=1;i<nodes.size();i++)
	{	if(coeff[i]<0) cout << "+";
		if(coeff[i]==1.)
			cout << "-";
		else if(coeff[i]!=-1.)
			cout << -coeff[i] << "*";
		cout << "nd[" << nodes[i] << "]";
	}
	if(qvalue>=0.) cout << "+";
	cout << qvalue << endl;
}
	

	
