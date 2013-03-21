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
    dir = ConvertToDirectionBits(dof);      // change to x,y,z bits
    angle1 = 0.;
    angle2 = 0.;
	
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
    dir = ConvertToDirectionBits(dof);      // change to x,y,z bits
    
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

// convert input dof to dir as bits x=1, y=2, z=4
int NodalVelBC::ConvertToDirectionBits(int dof)
{
    switch(dof)
    {   case X_DIRECTION:
        case Y_DIRECTION:
            return dof;
        case Z_DIRECTION_INPUT:
            return Z_DIRECTION;
        case XY_SKEWED_INPUT:
            return XY_SKEWED_DIRECTION;
        case XZ_SKEWED_INPUT:
            return XZ_SKEWED_DIRECTION;
        case YZ_SKEWED_INPUT:
            return YZ_SKEWED_DIRECTION;
        case XYZ_SKEWED_INPUT:
            return XYZ_SKEWED_DIRECTION;
        default:
            return 0;
    }
}

// convert dir bits to input dof
int NodalVelBC::ConvertToInputDof(void)
{
    int dof = dir&X_DIRECTION ? 1 : 0;
    if(dir&Y_DIRECTION)
    {   // will be 0 or 1
        dof = dof==0 ? 2 : 12 ;
    }
    if(dir&Z_DIRECTION)
    {   // will be 0, 1, 2, or 12
        if(dof==0)
            dof = 3;
        else if(dof==1)
            dof = 13;
        else if(dof==2)
            dof = 23;
        else
            dof = 123;
    }
    return dof;
}

#pragma mark NodalVelBC::Methods

// print it
BoundaryCondition *NodalVelBC::PrintBC(ostream &os)
{
    char nline[200];
	sprintf(nline,"%7d %3d %2d %15.7e %15.7e %7.2lf %7.2lf",nodeNum,ConvertToInputDof(),style,value,ftime,angle1,angle2);
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

// set to zero in x, y, or z velocity
NodalVelBC *NodalVelBC::ZeroVelBC(double mstime)
{	// set if has been activated
	int i = GetNodeNum(mstime);
	if(i>0) nd[i]->SetMomVel(dir);
    return (NodalVelBC *)GetNextObject();
}

// superpose x, y, or z velocity
NodalVelBC *NodalVelBC::AddVelBC(double mstime)
{	// set if has been activated
	int i = GetNodeNum(mstime);
	if(i>0)
	{	currentValue = BCValue(mstime);
		nd[i]->AddMomVel(dir,currentValue);
	}
    return (NodalVelBC *)GetNextObject();
}

// change to a ghost BC
NodalVelBC *NodalVelBC::SetGhostVelBC(double mstime)
{   // this will need to by BC property
    // distance in nd[] array to neighbors in direction dir, verified stays in grid
    int ghost = -1;         

    // set if has been activated
	int i = GetNodeNum(mstime);
	if(i>0 && nd[i]->NumberParticles())
	{	// see if neighbor in ghost direction fixes same dof
		if(nd[i+ghost]->fixedDirection&dir)
		{	// second node must by unfixed and have points
			int mirror = i+2*ghost;
			if(nd[mirror]->fixedDirection==0 && nd[mirror]->NumberParticles()>0)
			{	// found node to mirror
				//cout << "# node " << mirror << " vs. " << i ;
                // get CM mass from pk on node mirror
			}
		}
	}
	return (NodalVelBC *)GetNextObject();
}

// superpose x, y, or z velocity
NodalVelBC *NodalVelBC::InitFtot(double mstime)
{	// set if has been activated
	int i = GetNodeNum(mstime);
	if(i>0) nd[i]->SetFtot(dir,timestep);
	return (NodalVelBC *)GetNextObject();
}

// superpose x, y, or z velocity
NodalVelBC *NodalVelBC::AddFtot(double mstime)
{	// set if has been activated
	int i = GetNodeNum(mstime);
	if(i>0)
	{	// use currentValue set earlier in this step
		nd[i]->AddFtot(dir,timestep,currentValue);
	}
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
		nextBC = nextBC->ZeroVelBC(mstime);
    
    // Now add all velocities to nodes with velocity BCs
    nextBC=firstVelocityBC;
    while(nextBC!=NULL)
		nextBC = nextBC->AddVelBC(mstime);
	
	// check for ghosts
    /*
    nextBC=firstVelocityBC;
    while(nextBC!=NULL)
		nextBC = nextBC->SetGhostVelBC(mstime);
    */
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
		nextBC = nextBC->InitFtot(mstime);
    
    // Now add each superposed velocity BC at incremented time
    nextBC=firstVelocityBC;
    while(nextBC!=NULL)
		nextBC = nextBC->AddFtot(mstime);
}

