/********************************************************************************
	NodalVelGradBC.cpp
	nairn-mpm-fea

	Created by John Nairn on Feb 27, 2019.
	Copyright (c) 2019 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Boundary_Conditions/NodalVelGradBC.hpp"
#include "Read_XML/Expression.hpp"
#include "Nodes/NodalPoint.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Exceptions/CommonException.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "System/UnitsController.hpp"

#pragma mark NodalVelGradBC::Constructors and Destructors

// MPM Constructors
// Warning: dof must be 1,2,3 only
NodalVelGradBC::NodalVelGradBC(int num,int dof,int setStyle,double argTime,double gradDepth)
		: NodalVelBC(num,dof,setStyle,0.,argTime,0.,0.)
{
	// force to be +1 or -1
	depth = fabs(gradDepth);
	side = gradDepth>=0. ? 1 : -1;
	firstNode = -1;
	gradFunction = NULL;
	dispFunction = NULL;
	extractGrad = false;
}

// Destructor (and it is virtual)
NodalVelGradBC::~NodalVelGradBC()
{	if(gradFunction!=NULL) delete gradFunction;
	if(dispFunction!=NULL) delete dispFunction;
}

// Velocity gradient BCs current only allow directions 1 to 3
// throws CommonException()
int NodalVelGradBC::ConvertToDirectionBits(int dof)
{
	int maxDof = fmobj->IsThreeD() ? 3 : 2 ;
	if(dof<1 || dof>maxDof)
	{	throw CommonException("Invalid velocity gradietn input direction (should be 1,2, or 3).",
							  "NodalVelGradBC::ConvertToDirectionBits");
	}
	return NodalVelBC::ConvertToDirectionBits(dof);
}

#pragma mark NodalVelBC::Methods

// print it
BoundaryCondition *NodalVelGradBC::PrintBC(ostream &os)
{
	// print full details, even though not used here
	char nline[200];
    size_t nlsize=200;
    snprintf(nline,nlsize,"%7d %3d %2d %15.7e %15.7e %7.2lf %7.2lf",nodeNum,ConvertToInputDof(),style,
			GetBCValueOut(),GetBCFirstTimeOut(),angle1,angle2);
	os << nline;
	
	os << " ";
	if(function!=NULL)
		os << function->GetString();
	else
		os << "0";
	os << " : " << (double)side*depth << " : ";
	if(gradFunction!=NULL)
		os << gradFunction->GetString();
	else if(extractGrad)
		os << "mirror";
	else
		os << "0";
	os << " : ";
	if(dispFunction!=NULL)
		os << dispFunction->GetString();
	else
		os << "0";
	os << endl;
	return (BoundaryCondition *)GetNextObject();
}

// Get BC values for this time step
NodalVelGradBC *NodalVelGradBC::GetCurrentBCValue(double bctime)
{	// set if has been activated
	int i = GetNodeNum(bctime);
	if(i<=0)
	{	currentValue = gradValue = position = 0.;
		return (NodalVelGradBC *)GetNextObject();
	}
	
	// get velocity
	currentValue = BCValue(bctime);
	
	// get position and starting node
	BCPosition(bctime,i);
	
	return (NodalVelGradBC *)GetNextObject();
}

// set to zero in x, y, or z velocity
NodalVelGradBC *NodalVelGradBC::ZeroVelocityBC(double bctime,int passType)
{	// set if has been activated
	int i = GetNodeNum(bctime);
	if(i<=0) return (NodalVelGradBC *)GetNextObject();

	// get the gradient
	// 1. Here instead of in GetCurrentBCValue() because extracted gradient depends on passType
	// 2. Here to be before extrapolated values on firstNode are changed by zeroing
	BCGradValue(bctime,i,passType);
	
	// currently getting reaction force from all nodes on this BC
	if(passType==GRID_FORCES_CALL) ZeroVector(&freaction);
	
	// Look from firstNode to 2 past wall (or first inactive)
	int bcNode = firstNode;
	int numPast = 0;
	double dist;
	while(true)
	{	// get distance to the node
		if(dir&X_DIRECTION)
			dist = nd[bcNode]->x - position;
		else if(dir&Y_DIRECTION)
			dist = nd[bcNode]->y - position;
		else
			dist = nd[bcNode]->z - position;
		if(nd[bcNode]->NodeHasParticles())
		{	// zero the node
			nd[bcNode]->ZeroVelocityBC(&norm,passType,timestep,&freaction);
			
			// exit if >2 active nodes beyond the wall (may be wrongly placed)
			if((double)side*dist>0.)
			{	numPast++;
				if(numPast>2) break;
			}
		}
		else if((double)side*dist>0.)
		{	// exit on first inactive node beyond the wall
			break;
		}
		
		// next node, but nodes beyond the limit are inactive or don't exist
		bcNode += nodeStep;
		if(side<0)
		{	if(bcNode<lastNode) break;
		}
		else
		{	if(bcNode>lastNode) break;
		}
	}
	
	return (NodalVelGradBC *)GetNextObject();
}

// set to zero in x, y, or z velocity
NodalVelGradBC *NodalVelGradBC::AddVelocityBC(double bctime,int passType)
{	// set if has been activated
	int i = GetNodeNum(bctime);
	if(i<=0)
	{	currentValue = gradValue = position = 0.;
		return (NodalVelGradBC *)GetNextObject();
	}
	
	// Look from firstNode to 2 past wall (or first inactive)
	int bcNode = firstNode;
	int numPast = 0;
	double dist;
	while(true)
	{	// get distance to the node
		if(dir&X_DIRECTION)
			dist = nd[bcNode]->x - position;
		else if(dir&Y_DIRECTION)
			dist = nd[bcNode]->y - position;
		else
			dist = nd[bcNode]->z - position;
		if(nd[bcNode]->NodeHasParticles())
		{	// velocity at node's position
			double viBC = currentValue + gradValue*dist;
			
			nd[bcNode]->AddVelocityBC(&norm,viBC,passType,timestep,&freaction);

			// exit if >2 active nodes beyond the wall (may be wrongly placed)
			if((double)side*dist>0.)
			{	numPast++;
				if(numPast>2) break;
			}
		}
		else if((double)side*dist>0.)
		{	// exit on first inactive node beyond the wall
			break;
		}
		
		// next node, but nodes beyond the limit are inactive or don't exist
		bcNode += nodeStep;
		if(side<0)
		{	if(bcNode<lastNode) break;
		}
		else
		{	if(bcNode>lastNode) break;
		}
	}
	
	return (NodalVelGradBC *)GetNextObject();
}

#pragma mark NodalVelGradBC: Accessors

// set gradient function (if forGradient is true) or displacement function (if false)
// throws std::bad_alloc, SAXException()
void NodalVelGradBC::SetGradFunction(char *bcFunction,bool forGradient)
{
	if(bcFunction==NULL)
	{	ThrowSAXException("Velocity gradient function of time and position is missing");
		return;
	}
	if(strlen(bcFunction)==0)
	{	ThrowSAXException("Velocity gradient function function of time and position is missing");
		return;
	}
	
	if(forGradient)
	{	if(gradFunction!=NULL && extractGrad)
		{	ThrowSAXException("Duplicate velocity gradient function of time and position");
			return;
		}
		
		if(strcmp(bcFunction,"mirror")==0)
			extractGrad = true;
		else
			gradFunction =  Expression::CreateExpression(bcFunction,"Velocity gradient function not valid");
	}
	else
	{	if(dispFunction!=NULL)
		{	ThrowSAXException("Duplicate velocity displacement function of time and position");
			return;
		}
		
		dispFunction =  Expression::CreateExpression(bcFunction,"Velocity displacement function not valid");
	}
}

// find current position of the boundary condition
// never call unless beyond the start time
void NodalVelGradBC::BCPosition(double stepTime,int ni)
{
	// start on the node
	if(dir&X_DIRECTION)
		position = nd[ni]->x;
	else if(dir&Y_DIRECTION)
		position = nd[ni]->y;
	else
		position = nd[ni]->z;

	// done if no function
	if(dispFunction==NULL)
	{	// for fixed wall, find first node just once
		if(firstNode<0)
		{	nodeStep=1;
			firstNode = mpmgrid.FindShiftedNodeFromNode(ni,position,dir,side,nodeStep,depth);
			lastNode = mpmgrid.FindLastNode(firstNode,dir,side);
		}
		return;
	}
	
	// evaluate function (note that position here is of the initial node)
	// Legacy units: function initial node positions in mm, time in ms, function should return mm
	double bsStart = GetBCFirstTime();
    // (see Expression vmap)
	double vars[7];
	vars[0] = 6.5;
	vars[1] = UnitsController::Scaling(1.e3)*(stepTime-bsStart);		// t-tstart
	GetPositionVars(vars);
	position += dispFunction->EvaluateFunction(vars);
	
	// find involved nodes
	nodeStep=1;
	firstNode = mpmgrid.FindShiftedNodeFromNode(ni,position,dir,side,nodeStep,depth);
	lastNode = mpmgrid.FindLastNode(firstNode,dir,side);
}

// evaluate function to get boundary condition value
// position must aready be set
void NodalVelGradBC::BCGradValue(double stepTime,int ni,int bcType)
{
	// if a provided function, evaluate it
	if(gradFunction!=NULL)
	{	// function is current position along dir (in mm) and other positions for initial node
		// time is in ms. Evaulated results should be in 1/sec
		double bsStart = GetBCFirstTime();
        // (see Expression vmap)
		double vars[7];
		vars[0] = 6.5;
		vars[1] = UnitsController::Scaling(1.e3)*(stepTime-bsStart);		// t-tstart
		GetPositionVars(vars);
		// Change moving direction to current position (see Expression vmap)
		if(dir&X_DIRECTION)
			vars[2] = position;
		else if(dir&Y_DIRECTION)
			vars[3] = position;
		else
			vars[4] = position;
		gradValue = gradFunction->EvaluateFunction(vars);
		return;
	}
	else if(!extractGrad)
	{	// no function and no mirroring - i.e., use zero gradient
		gradValue = 0.;
		return;
	}
	
	// extract gradient from nodes firstNode-nodeStep and firstNode
	gradValue = 0.;
	
	// get c-o-m velocities, nodes firstNode-nodeStep to wall (exit if any empty)
	bool useVelocity = bcType==UPDATE_GRID_STRAINS_CALL ? true : false ;
	
	// least squares fit for node firstNode-nodeStep and firstNode to wall
	// I tried fitting all points before the wall, but seemed slightly worse
	int mnode = firstNode-nodeStep;
	Vector vi;
	double dvdx = 0., dxdx = 0., dx,dv;
	for(int nn=0;nn<2;nn++)
	{	// use zero gradient if any target nodes are inactive
		if(!nd[mnode]->GetCenterOfMassVelocity(&vi,useVelocity)) return;
		
		// find difference from current position value
		if(dir&X_DIRECTION)
		{	dx = nd[mnode]->x - position;
			dv = vi.x - currentValue;
		}
		else if(dir&Y_DIRECTION)
		{	dx = nd[mnode]->y - position;
			dv = vi.y - currentValue;
		}
		else
		{	dx = nd[mnode]->z - position;
			dv = vi.z - currentValue;
		}
		
		// add this point
		dxdx += dx*dx;
		dvdx += dv*dx;
		
		// next node
		mnode += nodeStep;
	}
	
	// least squares line for mirrored internal nodes and through currentValue at position
	gradValue = dvdx/dxdx;
}
