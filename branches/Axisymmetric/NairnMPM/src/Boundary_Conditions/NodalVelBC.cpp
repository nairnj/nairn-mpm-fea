/********************************************************************************
    NodalVelBC.cpp
    NairnMPM
    
    Created by John Nairn on Wed Jan 31 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Nodes/NodalPoint.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Exceptions/CommonException.hpp"

// Nodal velocity BC globals
NodalVelBC *firstVelocityBC=NULL;
NodalVelBC *lastVelocityBC=NULL;
NodalVelBC *firstRigidVelocityBC=NULL;
NodalVelBC *reuseRigidVelocityBC=NULL;

#pragma mark NodalVelBC::Constructors and Destructors

// MPM Constructors
NodalVelBC::NodalVelBC(int num,int dof,int setStyle,double velocity,double argTime)
		: BoundaryCondition(setStyle,velocity,argTime)
{
    nodeNum=num;
    dir = dof==3 ? Z_DIRECTION : dof ;		// change 3 to Z_DIRECTION (4) bit location
	
    // old dir==0 was skwed condition, now do by setting two velocities, thus never 0 here
    nd[nodeNum]->SetFixedDirection(dir);		// x, y, or z (1,2,4) directions
	
	pk=NULL;
}

// Destructor (and it is virtual)
NodalVelBC::~NodalVelBC()
{	if(pk!=NULL) delete pk;
}

// Reuse Rigid properties
BoundaryCondition *NodalVelBC::SetRigidProperties(int num,int dof,int setStyle,double velocity)
{	// set dir and direction
    dir = dof==3 ? Z_DIRECTION : dof ;		// change 3 to Z_DIRECTION (4) bit location
    
    // old dir==0 was skwed condition, now do by setting two velocities, thus never 0 here
    nd[num]->SetFixedDirection(dir);		// x, y, or z (1,2,4) directions
	
	// finish in base class (nodenum set there)
	return BoundaryCondition::SetRigidProperties(num,dof,setStyle,velocity);
}

// just unset condition, because may want to reuse it, return next one to unset
BoundaryCondition *NodalVelBC::UnsetDirection(void)
{	nd[nodeNum]->UnsetFixedDirection(dir);		// x, y, or z (1,2,4) directions
	return (BoundaryCondition *)GetNextObject();
}

#pragma mark NodalVelBC::Methods

// print it
BoundaryCondition *NodalVelBC::PrintBC(ostream &os)
{
    char nline[200];
	int outdir = dir==Z_DIRECTION ? 3 : dir ;
	sprintf(nline,"%7d %2d %2d %15.7e %15.7e",nodeNum,outdir,style,value,ftime);
    os << nline;
	PrintFunction(os);
	return (BoundaryCondition *)GetNextObject();
}

// save nodal momentum and velocity (but only once)
NodalVelBC *NodalVelBC::CopyNodalVelocities(NodalPoint *nd)
{
	// create vector to hold options
#ifdef _BC_CRACK_SIDE_ONLY_
	if(pk==NULL)
	{	pk=(Vector *)malloc(sizeof(Vector)*maxMaterialFields);
		if(pk==NULL) throw CommonException("Memory error allocating vectors to copy boundary condition nodal values.",
										   "NodalVelBC::CopyNodalVelocities");
	}
	nd->cvf[0]->CopyFieldMomenta(pk,0);
#else
	if(pk==NULL)
	{	pk=(Vector *)malloc(sizeof(Vector)*maxMaterialFields*maxCrackFields);
		if(pk==NULL) throw CommonException("Memory error allocating vectors to copy boundary condition nodal values.",
										   "NodalVelBC::CopyNodalVelocities");
	}
	int i;
	int offset=0;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(nd->cvf[i]))
			offset=nd->cvf[i]->CopyFieldMomenta(pk,offset);
	}
#endif
    return (NodalVelBC *)GetNextObject();
}

// paste nodal momentum and velocity (but only once)
NodalVelBC *NodalVelBC::PasteNodalVelocities(NodalPoint *nd)
{
#ifdef _BC_CRACK_SIDE_ONLY_
	nd->cvf[0]->PasteFieldMomenta(pk,0);
#else
	int i;
	int offset=0;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveField(nd->cvf[i]))
			offset=nd->cvf[i]->PasteFieldMomenta(pk,offset);
	}
#endif
    return (NodalVelBC *)GetNextObject();
}

#pragma mark NodelVelBC::Class Methods

/*******************************************************************
    Impose specified momenta at selected nodes.
    The imposed momenta are needed before any strain update.
	Called in Tasks 1 and 6
    
    Note: makeCopy is TRUE for Task 1 and FALSE for Task 6
        In Task 6, use BC at mtime+timestep
    
    When makeCopy is true (Task 1), it stores
        copy of the original nodal momenta that correspond to 
        initial particle extrapolation. After calculating total
        force by extrapolation, these values are pasted back
        and the forces are changed such that the momentum
        update will result in nodal momenta specified by
        the boundary conditions.
    
    This approach is needed such that the nodal forces (which
        are essentially accelerations) will have the right
        values to correctly update particle velocities. In other
		words, the accelerations are adjusted to cause the
		no-boundary-condition momenta to update to the
		specifed momenta. Without this fix the particle positions
		would be correct (velocities are OK), but the particle
		velocities would be wrong (accelerations not consistent).
*/
void NodalVelBC::GridMomentumConditions(int makeCopy)
{
    int i;
    NodalVelBC *nextBC;
	double mstime;
    
	// convert time to ms, use time and beginning or end of time step
    // On first pass (when true), copy nodal momenta before anything changed
    //	(may make multiple copies, but that is OK)
	if(makeCopy)
	{	mstime=1000.*mtime;
    	nextBC=firstVelocityBC;
        while(nextBC!=NULL)
		{	i=nextBC->GetNodeNum();
			nextBC=nextBC->CopyNodalVelocities(nd[i]);
        }
    }
	else
		mstime=1000.*(mtime+timestep);
		
    // Now zero nodes with velocity set by BC
    nextBC=firstVelocityBC;
    while(nextBC!=NULL)
    {	// x, y, or z velocity will be set
        if((i=nextBC->GetNodeNum(mstime)))
            nd[i]->SetMomVel(nextBC->dir);
        nextBC=(NodalVelBC *)nextBC->GetNextObject();
    }
    
    // Now add all velocities to nodes with velocity BCs
    nextBC=firstVelocityBC;
    while(nextBC!=NULL)
	{	// x, y, or z velocity is incremented
		if((i=nextBC->GetNodeNum(mstime)))
            nd[i]->AddMomVel(nextBC->dir,nextBC->BCValue(mstime));
        nextBC=(NodalVelBC *)nextBC->GetNextObject();
    }
}

/**********************************************************
    Set forces on nodes with set velocity such that
    momentum after the update will match the BC
    velocity (and momentum)
*/
void NodalVelBC::ConsistentGridForces(void)
{
    int i;
    NodalVelBC *nextBC=firstVelocityBC;
    
	// advanced time in msec
	double mstime=1000.*(mtime+timestep);
	
    // First restore initial nodal values
    while(nextBC!=NULL)
	{	i=nextBC->GetNodeNum();
		nextBC=nextBC->PasteNodalVelocities(nd[i]);
    }
    
    // Second set force to -p(interpolated)/timestep
    nextBC=firstVelocityBC;
    while(nextBC!=NULL)
	{	// x velocity will be set
		if((i=nextBC->GetNodeNum(mstime)))
            nd[i]->SetFtot(nextBC->dir,timestep);
        nextBC=(NodalVelBC *)nextBC->GetNextObject();
    }
    
    // Now add each superposed velocity BC at incremented time
    nextBC=firstVelocityBC;
    while(nextBC!=NULL)
	{	// x velocity will be set
		if((i=nextBC->GetNodeNum(mstime)))
            nd[i]->AddFtot(nextBC->dir,timestep,nextBC->BCValue(mstime));
        nextBC=(NodalVelBC *)nextBC->GetNextObject();
    }
}

