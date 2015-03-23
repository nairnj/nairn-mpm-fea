/********************************************************************************
    NodalVelBC.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 31 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Nodes/NodalPoint.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Exceptions/CommonException.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"

// Nodal velocity BC globals
NodalVelBC *firstVelocityBC=NULL;
NodalVelBC *lastVelocityBC=NULL;
NodalVelBC *firstRigidVelocityBC=NULL;
NodalVelBC *reuseRigidVelocityBC=NULL;

#pragma mark NodalVelBC::Constructors and Destructors

// MPM Constructors
// Warning: dof must be Z_DIRECTION_INPUT (=3) to get Z axis velocity and not Z_DIRECTION (=4)
//          input numbers are 1,2,3,12,13,23,123
NodalVelBC::NodalVelBC(int num,int dof,int setStyle,double velocity,double argTime,double ang1,double ang2)
		: BoundaryCondition(setStyle,velocity,argTime)
{
    nodeNum = num;
	reflectedNode = -1;
    mirrorSpacing = 0;                      // +1 shift higher nodes, -1 shift lower, 0 no mirroring
    angle1 = ang1;
    angle2 = ang2;
    dir = ConvertToDirectionBits(dof);      // change input settings to x,y,z bits
    SetNormalVector();						// get normal vector depending on dir and angles
    ZeroVector(&freaction);
	
    // old dir==0 was skwed condition, now do by setting two velocities, thus never 0 here
    nd[nodeNum]->SetFixedDirection(dir);		// x, y, or z (1,2,4) directions
	
	pk=NULL;
}

// Destructor (and it is virtual)
NodalVelBC::~NodalVelBC()
{	if(pk!=NULL) delete pk;
}

// Reuse Rigid properties
// Warning: dof must be BC type and not input type
BoundaryCondition *NodalVelBC::SetRigidProperties(int num,int dof,int setStyle,double velocity)
{	// set dir and direction (cannot yet be skewed directions)
    angle1 = 0.;
    angle2 = 0.;
    dir = dof;                          // internal calls already have correct bit settings
    SetNormalVector();                  // get normal vector depending on dir and angles
    
    // old dir==0 was skewed condition, now do by setting two velocities, thus never 0 here
    nd[num]->SetFixedDirection(dir);		// x, y, or z (1,2,4) directions
    
    // not reflected yet
    reflectedNode = -1;
    mirrorSpacing = 0;                      // +1 shift higher nodes, -1 shift lower, 0 no mirroring
	
	// finish in base class (nodenum set there)
	return BoundaryCondition::SetRigidProperties(num,dof,setStyle,velocity);
}

// just unset condition, because may want to reuse it, return next one to unset
BoundaryCondition *NodalVelBC::UnsetDirection(void)
{	nd[nodeNum]->UnsetFixedDirection(dir);		// x, y, or z (1,2,4) directions
	return (BoundaryCondition *)GetNextObject();
}

// convert input dof to dir as bits x=1, y=2, z=4
// input options as 1,2,3,12,13,23,123 and changes to bitwise settings
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
            throw CommonException("Invalid velocity input direction (should be 1,2,3,12,13,23, or 123).",
                                  "NodalVelBC::ConvertToDirectionBits");
            return 0;
    }
}

// initialize mormal vector depending on direction bits
// dir is 1 to 7 (3 axis bits can be set)
void NodalVelBC::SetNormalVector(void)
{
    switch(dir)
    {   case X_DIRECTION:
			norm = MakeVector(1.,0.,0.);
            break;
        case Y_DIRECTION:
			norm = MakeVector(0.,1.,0.);
            break;
        case Z_DIRECTION:
			norm = MakeVector(0.,0.,1.);
            break;
        case XY_SKEWED_DIRECTION:
			norm = MakeVector(cos(PI_CONSTANT*angle1/180.),-sin(PI_CONSTANT*angle1/180.),0.);
            break;
        case XZ_SKEWED_DIRECTION:
			norm = MakeVector(cos(PI_CONSTANT*angle1/180.),0,-sin(PI_CONSTANT*angle1/180.));
            break;
        case YZ_SKEWED_DIRECTION:
			norm = MakeVector(0.,cos(PI_CONSTANT*angle1/180.),-sin(PI_CONSTANT*angle1/180.));
            break;
        case XYZ_SKEWED_DIRECTION:
			norm = MakeVector(cos(PI_CONSTANT*angle2/180.)*sin(PI_CONSTANT*angle1/180.),
                              sin(PI_CONSTANT*angle2/180.)*sin(PI_CONSTANT*angle1/180.),
                              cos(PI_CONSTANT*angle1/180.));
            break;
        default:
            throw CommonException("Invalid direction bits (should be 0x000 to 0x111).",
                                  "NodalVelBC::SetNormalVector");
            break;
    }
 }


// convert dir bits to input dof
// only used with printing velocity BCs to output stream
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
	sprintf(nline,"%7d %3d %2d %15.7e %15.7e %7.2lf %7.2lf",nodeNum,ConvertToInputDof(),style,
							GetBCValueOut(),GetBCFirstTimeOut(),angle1,angle2);
    os << nline;
	PrintFunction(os);
	return (BoundaryCondition *)GetNextObject();
}

// save nodal momentum and velocity (but only once)
NodalVelBC *NodalVelBC::CopyNodalVelocities(NodalPoint *nd)
{
	// create vector to hold options
	if(pk==NULL)
	{	pk=(Vector *)malloc(sizeof(Vector)*maxMaterialFields*maxCrackFields);
		if(pk==NULL) throw CommonException("Memory error allocating vectors to copy boundary condition nodal values.",
										   "NodalVelBC::CopyNodalVelocities");
	}
	int i;
	int offset=0;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveNonrigidField(nd->cvf[i]))
			offset=nd->cvf[i]->CopyFieldMomenta(pk,offset);
	}
    
    return (NodalVelBC *)GetNextObject();
}

// paste nodal momentum and velocity (but only once)
NodalVelBC *NodalVelBC::PasteNodalVelocities(NodalPoint *nd)
{
	int i;
	int offset=0;
	for(i=0;i<maxCrackFields;i++)
	{	if(CrackVelocityField::ActiveNonrigidField(nd->cvf[i]))
			offset=nd->cvf[i]->PasteFieldMomenta(pk,offset);
	}
    return (NodalVelBC *)GetNextObject();
}

// set to zero in x, y, or z velocity
NodalVelBC *NodalVelBC::ZeroVelBC(double bctime)
{	// set if has been activated
	int i = GetNodeNum(bctime);
	if(i>0) nd[i]->SetMomVel(&norm);
    return (NodalVelBC *)GetNextObject();
}

// superpose x, y, or z velocity
NodalVelBC *NodalVelBC::AddVelBC(double bctime)
{	// set if has been activated
	int i = GetNodeNum(bctime);
	if(i>0)
    {   currentValue = BCValue(bctime);
		if(reflectedNode<0)
		{	// scalar value in norm direction
			nd[i]->AddMomVel(&norm,currentValue);
		}
		else
		{	// reflect one component at a symmetry plane
            nd[i]->AddMomVel(&norm,2.*currentValue);
			nd[i]->ReflectMomVel(&norm,nd[reflectedNode]);
		}
	}
    return (NodalVelBC *)GetNextObject();
}

// For rigid BC, see if it has a refleceted node
NodalVelBC *NodalVelBC::SetMirroredVelBC(double bctime)
{
    // exit if not doing mirroring
    if(mirrorSpacing==0) return (NodalVelBC *)GetNextObject();
    
    // set if has been activated
	int i = GetNodeNum(bctime);
    if(i>0)
    {   // look at neighbor in mirror direction if node i has particles and neighbor is in the grid
        int neighbor = i+mirrorSpacing;
        if(nd[i]->NodeHasNonrigidParticles() && neighbor>0 && neighbor<=nnodes)
        {	// see if neighbor fixes same dof
            if(nd[neighbor]->fixedDirection&dir)
            {	// the next node in mirror direction if in the grid
                int mirror = neighbor+mirrorSpacing;
                if(mirror>0 && mirror<=nnodes)
                {   // it should not fix the direection and should have particles
                    if((nd[mirror]->fixedDirection&dir)==0 && nd[mirror]->NodeHasNonrigidParticles()>0)
                    {   // maybe should verify neighbor is for same rigid material, but needs to check all BCs for their ID
                        
                        // found node to reflect
                        //cout << "#    node " << mirror << " from " << i << " through " << neighbor << endl;
                        reflectedNode = mirror;
                    }
                }
            }
        }
	}
	return (NodalVelBC *)GetNextObject();
}

// Initialize ftot to -(pk.norm/deltime) norm in each material velocity field
// freaction will be sum over all material velocity fields for this BC only
// freaction will be zero for second BC on same node with same norm
// total reaction on node us sum of freaction over all BCs
NodalVelBC *NodalVelBC::InitFtotDirection(double bctime)
{	// set if has been activated
	int i = GetNodeNum(bctime);
    ZeroVector(&freaction);
	if(i>0) nd[i]->SetFtotDirection(&norm,timestep,&freaction);
	return (NodalVelBC *)GetNextObject();
}

// superpose force in direction normal
// add force for this BC to freaction
NodalVelBC *NodalVelBC::SuperposeFtotDirection(double bctime)
{	// set if has been activated
	int i = GetNodeNum(bctime);
	if(i>0)
	{	if(reflectedNode<0)
		{	// use currentValue set earlier in this step
			nd[i]->AddFtotDirection(&norm,timestep,currentValue,&freaction);
		}
		else
		{	// use reflected velocity
            nd[i]->AddFtotDirection(&norm,timestep,2.*currentValue,&freaction);
			nd[i]->ReflectFtotDirection(&norm,timestep,nd[reflectedNode],&freaction);
		}
	}
	return (NodalVelBC *)GetNextObject();
}

// when getting total reaction force, add freaction to input vector
// if matchID==0 include it, otherwise include only in matchID equals bcID
NodalVelBC *NodalVelBC::AddReactionForce(Vector *totalReaction,int matchID)
{	if(bcID==matchID || matchID==0)
		AddVector(totalReaction,&freaction);
	return (NodalVelBC *)GetNextObject();
}

// On symmetry plane set velocity to minus component of the reflected node
void NodalVelBC::SetReflectedNode(int mirrored) { reflectedNode = mirrored; }

// On rigid particle set spacing to mirrored node
// Needs structured grid, but does not appear to need equal element sizes
void NodalVelBC::SetMirrorSpacing(int mirrored)
{
	mirrorSpacing = 0;
	if(mirrored==0) return;
	if(!mpmgrid.IsStructuredGrid()) return;
	
	switch(dir)
    {   case X_DIRECTION:
			mirrorSpacing = mirrored<0 ? mpmgrid.xplane : -mpmgrid.xplane ;
            break;
        case Y_DIRECTION:
			mirrorSpacing = mirrored<0 ? mpmgrid.yplane : -mpmgrid.yplane ;
            break;
        case Z_DIRECTION:
			mirrorSpacing = mirrored<0 ? mpmgrid.zplane : -mpmgrid.zplane ;
            break;
		default:
			break;
	}

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
    
	// convert time to ms, use time at beginning or end of time step
    // On first pass (when true), copy nodal momenta before anything changed
    //	(may make multiple copies, but that is OK)
	if(makeCopy)
	{	nextBC=firstVelocityBC;
        while(nextBC!=NULL)
		{	i=nextBC->GetNodeNum(mtime);
			if(i>0) nextBC->CopyNodalVelocities(nd[i]);
			nextBC = (NodalVelBC *)nextBC->GetNextObject();
        }
    }
	
    // Now zero nodes with velocity set by BC
    nextBC=firstVelocityBC;
    while(nextBC!=NULL)
		nextBC = nextBC->ZeroVelBC(mtime);
    
    // Now add all velocities to nodes with velocity BCs
    nextBC=firstVelocityBC;
    while(nextBC!=NULL)
		nextBC = nextBC->AddVelBC(mtime);
}

/**********************************************************
    Set forces on nodes with set velocity such that
    momentum after the update will match the BC
    velocity (and momentum)
*/
void NodalVelBC::ConsistentGridForces(void)
{
    int i;
    NodalVelBC *nextBC = firstVelocityBC;
    
    // First restore initial nodal values
    while(nextBC!=NULL)
	{	i = nextBC->GetNodeNum(mtime);
		if(i>0) nextBC->PasteNodalVelocities(nd[i]);
		nextBC = (NodalVelBC *)nextBC->GetNextObject();
    }
	
    // Second set force to -p(interpolated)/timestep
    nextBC = firstVelocityBC;
    while(nextBC!=NULL)
		nextBC = nextBC->InitFtotDirection(mtime);
    
	// Now add each superposed velocity BC at incremented time
    nextBC=firstVelocityBC;
    while(nextBC!=NULL)
		nextBC = nextBC->SuperposeFtotDirection(mtime);
	
}

/**********************************************************
	Sum all reaction forces for all velocity BCs
	If ID is not zero, only includes those with matching ID
	Result is in micro Newtons
*/
Vector NodalVelBC::TotalReactionForce(int matchID)
{
	Vector reactionTotal;
	ZeroVector(&reactionTotal);
    NodalVelBC *nextBC=firstVelocityBC;
    
    // First restore initial nodal values
    while(nextBC!=NULL)
		nextBC=nextBC->AddReactionForce(&reactionTotal,matchID);
	
	return reactionTotal;
}
	


