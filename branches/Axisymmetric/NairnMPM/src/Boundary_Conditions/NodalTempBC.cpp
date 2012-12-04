/********************************************************************************
    NodalTempBC.cpp
    NairnMPM
    
    Created by John Nairn on Oct 18 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Boundary_Conditions/NodalTempBC.hpp"
#include "Nodes/NodalPoint.hpp"

// Nodal temperature BC global
NodalTempBC *firstTempBC=NULL;
NodalTempBC *lastTempBC=NULL;
NodalTempBC *firstRigidTempBC=NULL;
NodalTempBC *reuseRigidTempBC=NULL;

/*******************************************************************
	NodalTempBC: Constructors and Destructors
*******************************************************************/

NodalTempBC::NodalTempBC(int num,int setStyle,double temperature,double argTime)
		: BoundaryCondition(setStyle,temperature,argTime)
{
    nodeNum=num;
	nd[nodeNum]->SetFixedDirection(TEMP_DIRECTION);
}

// Reuse Rigid BC
BoundaryCondition *NodalTempBC::SetRigidProperties(int num,int dof,int setStyle,double temperature)
{	// set direction
	nd[num]->SetFixedDirection(TEMP_DIRECTION);
	// finish in base class
	return BoundaryCondition::SetRigidProperties(num,dof,setStyle,temperature);
}

// just unset condition, because may want to reuse it, return next one to unset
BoundaryCondition *NodalTempBC::UnsetDirection(void)
{	nd[nodeNum]->UnsetFixedDirection(TEMP_DIRECTION);
	return (BoundaryCondition *)GetNextObject();
}

/*******************************************************************
	NodalTempBC: Methods
*******************************************************************/

// print it
BoundaryCondition *NodalTempBC::PrintBC(ostream &os)
{
    char nline[200];
	sprintf(nline,"%7d %2d %15.7e %15.7e",nodeNum,style,value,ftime);
    os << nline;
	PrintFunction(os);
	return (BoundaryCondition *)GetNextObject();
}

// save nodal temperature and zero it
NodalTempBC *NodalTempBC::CopyNodalTemperature(NodalPoint *nd)
{
	temperatureNoBC=nd->gTemperature;
	return (NodalTempBC *)GetNextObject();
}

// restore nodal temperature and get initial force to cancel no-BC result
NodalTempBC *NodalTempBC::PasteNodalTemperature(NodalPoint *nd)
{
	nd->gTemperature=temperatureNoBC;
	return (NodalTempBC *)GetNextObject();
}
