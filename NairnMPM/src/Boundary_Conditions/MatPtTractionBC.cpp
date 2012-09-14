/********************************************************************************
	MatPtTractionBC.cpp
	NairnMPM

	Created by John Nairn on 9/13/12.
	Copyright (c) 2012 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Boundary_Conditions/MatPtTractionBC.hpp"
#include "MPM_Classes/MPMBase.hpp"

// global
MatPtTractionBC *firstTractionPt=NULL;

#pragma mark MatPtTractionBC::Constructors and Destructors

// Constructors
MatPtTractionBC::MatPtTractionBC(int num,int dof,int edge,int sty)
							: MatPtLoadBC(num,dof,sty)
{
	face = edge;
}

#pragma mark MatPtTractionBC: Methods

// print to output
BoundaryCondition *MatPtTractionBC::PrintBC(ostream &os)
{
    char nline[200];
    
    sprintf(nline,"%7d %2d   %2d  %2d %15.7e %15.7e",ptNum,direction,face,style,value,ftime);
    os << nline;
	PrintFunction(os);
	
    return (BoundaryCondition *)GetNextObject();
}

// increment external load on a particle
// input is analysis time in seconds
MatPtTractionBC *MatPtTractionBC::AddMPTraction(double bctime)
{
	double mstime=1000.*bctime;
	MPMBase *mpmptr = mpm[ptNum-1];
	double tmag = BCValue(mstime);
	
	// get corners and direction from material point
	int cElem[4];
	Vector corners[4],tnorm;
	mpmptr->GetTractionInfo(face,direction,cElem,corners,&tnorm);
	
	// loop over corner finding all nodes and add to fext
		
    return (MatPtTractionBC *)GetNextObject();
}

#pragma mark MatPtTractionBC: Class Methods

// Calculate traction forces applied to particles and add to nodal fext
void MatPtTractionBC::SetParticleSurfaceTractions(double stepTime)
{
    MatPtTractionBC *nextLoad=firstTractionPt;
    while(nextLoad!=NULL)
    	nextLoad=nextLoad->AddMPTraction(stepTime);
}

