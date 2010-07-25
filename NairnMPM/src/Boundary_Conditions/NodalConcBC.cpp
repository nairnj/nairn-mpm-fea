/********************************************************************************
    NodalConcBC.cpp
    NairnMPM
    
    Created by John Nairn on Thu Apr 1 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Boundary_Conditions/NodalConcBC.hpp"
#include "Nodes/NodalPoint.hpp"

// Nodal concentration BC globals
NodalConcBC *firstConcBC=NULL;
NodalConcBC *lastConcBC=NULL;
NodalConcBC *firstRigidConcBC=NULL;
NodalConcBC *reuseRigidConcBC=NULL;

/*******************************************************************
	NodalConcBC: Constructors and Destructors
*******************************************************************/

NodalConcBC::NodalConcBC(long num,int setStyle,double concentration,double argTime)
		: BoundaryCondition(setStyle,concentration,argTime)
{
    nodeNum=num;
	nd[nodeNum]->SetFixedDirection(CONC_DIRECTION);
}

// Reuse Rigid BC
BoundaryCondition *NodalConcBC::SetRigidProperties(long num,int dof,int setStyle,double concentration)
{	// set direction
	nd[num]->SetFixedDirection(CONC_DIRECTION);
	// finish in base class
	return BoundaryCondition::SetRigidProperties(num,dof,setStyle,concentration);
}

// just unset condition, because may want to reuse it, return next one to unset
BoundaryCondition *NodalConcBC::UnsetDirection(void)
{	nd[nodeNum]->UnsetFixedDirection(CONC_DIRECTION);
	return (BoundaryCondition *)GetNextObject();
}

/*******************************************************************
	NodalConcBC: Methods
*******************************************************************/

// print it
BoundaryCondition *NodalConcBC::PrintBC(ostream &os)
{
    char nline[200];
	sprintf(nline,"%5ld %2d %15.7e %15.7e",nodeNum,style,value,ftime);
    os << nline;
	PrintFunction(os);
	return (BoundaryCondition *)GetNextObject();
}

// save nodal concentration and zero it
NodalConcBC *NodalConcBC::CopyNodalConcentration(NodalPoint *nd)
{
	concentrationNoBC=nd->gConcentration;
	return (NodalConcBC *)GetNextObject();
}

// restore nodal concentration and get initial force to cancel no-BC result
NodalConcBC *NodalConcBC::PasteNodalConcentration(NodalPoint *nd)
{
	nd->gConcentration=concentrationNoBC;
	return (NodalConcBC *)GetNextObject();
}
