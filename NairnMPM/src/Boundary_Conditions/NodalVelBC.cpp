/********************************************************************************
    NodalVelBC.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 31 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/
#if defined ( _MSC_VER) || defined (__APPLE__) 
#include "stdafx.h"
#endif
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Nodes/NodalPoint.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Exceptions/CommonException.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Global_Quantities/BodyForce.hpp"
#include "NairnMPM_Class/UpdateStrainsFirstTask.hpp"
#include <cmath>
#include "tools/ExtractMPM.hpp" //PI

// Nodal velocity BC globals
NodalVelBC *firstVelocityBC=NULL;
NodalVelBC *lastVelocityBC=NULL;
NodalVelBC *firstRigidVelocityBC=NULL;
NodalVelBC *reuseRigidVelocityBC=NULL;

bool NodalVelBC::holdAllVelocityBCs = false;

#pragma mark NodalVelBC::Constructors and Destructors

// MPM Constructors
// Warning: dof must be Z_DIRECTION_INPUT (=3) to get Z axis velocity and not Z_DIRECTION (=4)
//          input numbers are 1,2,3,12,13,23,123
NodalVelBC::NodalVelBC(int num,int dof,int setStyle,double velocity,double argTime,double ang1,double ang2)
		: BoundaryCondition(setStyle,velocity,argTime)
{
    nodeNum = num;
	reflectedNode = -1;
	reflectRatio = 1.;
    mirrorSpacing = 0;                      // +1 shift higher nodes, -1 shift lower, 0 no mirroring
    angle1 = ang1;
    angle2 = ang2;
    dir = ConvertToDirectionBits(dof);      // change input settings to x,y,z bits
    SetNormalVector();						// get normal vector depending on dir and angles
    ZeroVector(&freaction);
	
	// old dir==0 was skewed condition, now do by setting two velocities, thus never 0 here
	nd[nodeNum]->SetFixedDirection(dir);		// three bits (1-8) for x, y, z or skewed in 2 or three directions
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
    // direct set in parent lass		// x, y, or z (1,2,4) directions
    
    // not reflected yet
    reflectedNode = -1;
    mirrorSpacing = 0;                      // +1 shift higher nodes, -1 shift lower, 0 no mirroring
	
	// finish in base class (nodenum set there)
	return BoundaryCondition::SetRigidProperties(num,dir,setStyle,velocity);
}

// get set direction
int NodalVelBC::GetSetDirection(void) const { return dir; }

// convert input dof to dir as bits x=1, y=2, z=4
// input options as 1,2,3,12,13,23,123 and changes to bitwise settings
// throws CommonException()
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
// throws CommonException()
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
            throw CommonException("Invalid direction bits (should be 0x001 to 0x111).",
                                  "NodalVelBC::SetNormalVector");
            break;
    }
 }

// return pointer to the BC's normal vector
Vector *NodalVelBC::GetNormalVector(void) { return &norm; }

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
    size_t nlsize=200;
	snprintf(nline,nlsize,"%7d %3d %2d %15.7e %15.7e %7.2lf %7.2lf",nodeNum,ConvertToInputDof(),style,
							GetBCValueOut(),GetBCFirstTimeOut(),angle1,angle2);
    os << nline;
	PrintFunction(os);
	return (BoundaryCondition *)GetNextObject();
}

// Get BC value for this time step
// superpose x, y, or z velocity
NodalVelBC *NodalVelBC::GetCurrentBCValue(double bctime)
{	// set if has been activated
	int i = GetNodeNum(bctime);
	currentValue = i>0 ? BCValue(bctime) : 0. ;
	return (NodalVelBC *)GetNextObject();
}

// Set to zero in x, y, or z velocity
// For grid forces subtraction ((Ftot.norm) + (pk.norm)/deltime)*norm from current force
// freaction will be sum over all material velocity fields for this BC only
// freaction will be zero for second BC on same node with same norm
// total reaction on node is sum of freaction over all BCs
NodalVelBC *NodalVelBC::ZeroVelocityBC(double bctime,int passType)
{	// set if has been activated
	int i = GetNodeNum(bctime);
	if(passType==GRID_FORCES_CALL) ZeroVector(&freaction);
	if(i>0) nd[i]->ZeroVelocityBC(&norm,passType,timestep,&freaction);
    return (NodalVelBC *)GetNextObject();
}

// superpose x, y, or z velocity
NodalVelBC *NodalVelBC::AddVelocityBC(double bctime,int passType)
{	// set if has been activated
	int i = GetNodeNum(bctime);
	if(i>0)
    {	//currentValue = BCValue(bctime);
		if(reflectedNode<0)
		{	// scalar value in norm direction
			nd[i]->AddVelocityBC(&norm,currentValue,passType,timestep,&freaction);
		}
		else
		{	// reflect one component at a symmetry plane
			nd[i]->ReflectVelocityBC(&norm,nd[reflectedNode],currentValue,reflectRatio,passType,timestep,&freaction);
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
    {   // look at neighbor in mirror direction but only if node i has particles and neighbor is in the grid
        int neighbor = i+mirrorSpacing;
        if(nd[i]->NodeHasNonrigidParticles() && neighbor>0 && neighbor<=nnodes)
        {	// see if neighbor fixes same dof
            if(nd[neighbor]->fixedDirection&dir)
            {	// the next node in mirror direction if in the grid
                int mirror = neighbor+mirrorSpacing;

                if(mirror>0 && mirror<=nnodes)
                {   // it should not fix the direection and should have particles
                    if((nd[mirror]->fixedDirection&dir)==0 && nd[mirror]->NodeHasNonrigidParticles())
                    {   // maybe should verify neighbor is for same rigid material, but needs to check all BCs for their ID
                        
                        // found node to reflect
                        reflectedNode = mirror;
						reflectRatio = mpmgrid.GetCellRatio(nd[neighbor],dir,mirrorSpacing);
                    }
                }
            }
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
void NodalVelBC::SetReflectedNode(int mirrored,double cellRatio)
{	reflectedNode = mirrored;
	reflectRatio = cellRatio;
}

// On rigid particle set spacing to mirrored node
// Needs structured grid, but does not need equal element sizes
// mirrorSpacing is node spacing in plane for this BCs node to be refected it needs:
//	a. mirrorSpacing away is BC with same fixed direction and 2*mirrorSpacing away is free node in the object
void NodalVelBC::SetMirrorSpacing(int mirrored)
{
	mirrorSpacing = 0;
	if(mirrored==0) return;
	
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
 Get BC value once per time step after mass and momentum
 update and after dynamic grid BCs are set up
 This avoid recalculation of currentValue, which only matters
	when it is function
*/
void NodalVelBC::GridVelocityBCValues(void)
{
	// skip if no BCs or ot development mode that omits this calculation
	if(firstVelocityBC==NULL) return;
	
	// Now zero nodes with velocity set by BC
	NodalVelBC *nextBC=firstVelocityBC;
	while(nextBC!=NULL)
		nextBC = nextBC->GetCurrentBCValue(mtime);
}

/*****************************************************************************
 Impose specified momenta, force, or velocity at selected nodes.
 The passTypes are
	MASS_MOMENTUM_CALL: after momentum extrapolation, use lumped mass
		matrix method to change pk
	GRID_FORCES_CALL: After force extrapolation, use lumped method
		to get grid forces corresponding to lumped momentum change
	UPDATE_MOMENTUM_CALL: after momentum update, use lumped mass matrix
		method to change pk and ftot. It is usually small change
		caused by contact or mirroring (might not be needed)
	UPDATE_STRAINS_LAST_CALL: prior to USL+ or USAVG+ use lumped mass
		matrix method to change pk
	UPDATE_GRID_STRAINS_CALL: For FMPM(k>1) impose velocity BCs in vk[0]
 Be sure to precheck simulation wants this calculation BEFORE calling
*/
void NodalVelBC::GridVelocityConditions(int passType)
{
	// skip in no velocity BCs
	if(firstVelocityBC==NULL) return;
	
	switch(passType)
	{	case UPDATE_GRID_STRAINS_CALL:
			// Only called for FMPM(k>1) after grid velocities are found by XPIC
			// BC options 0 and 1 impose BCs in these velocties now
			VelocityBCLoop(UPDATE_GRID_STRAINS_CALL);
			break;
		
		case MASS_MOMENTUM_CALL:
		{
#if ADJUST_COPIED_PK == 1
			// I think this should only be done after initial extrapolation
			// adjust for symmetry plane option
			NodalVelBC *nextBC=firstVelocityBC;
			while(nextBC!=NULL)
			{	int i=nextBC->GetNodeNum();
				if(nd[i]->fixedDirection&ANYSYMMETRYPLANE_DIRECTION)
				{	for(int j=0;j<maxCrackFields;j++)
					{	if(CrackVelocityField::ActiveField(nd[i]->cvf[j]))
						nd[i]->cvf[j]->AdjustForSymmetryBC(nd[i]);
					}
				}
				nextBC = (NodalVelBC *)nextBC->GetNextObject();
			}
#endif
			// no need when no strain update being done
			if(USFTask==NULL) return;
			
			// Do the loop
			VelocityBCLoop(passType);
			break;
		}
		case UPDATE_STRAINS_LAST_CALL:
		case GRID_FORCES_CALL:
			// Here for grid forces and update strains last
			VelocityBCLoop(passType);
			break;
			
		case UPDATE_MOMENTUM_CALL:
			// after update, momenta should be about correct
			// may be better to repeat if contact or mirrored BCs
			// Not sure what FMPM should do yet, so does nothing for now
			if(bodyFrc.GetXPICOrder()<=1)
			{	// FLIP and FMPM(1)
				VelocityBCLoop(passType);
			}
			break;
		
		default:
			break;
	}
}

/**********************************************************
	Loop of BCs zeroing and then adding BCs
 	Loop using specified passType
*/
void NodalVelBC::VelocityBCLoop(int passType)
{
	// Now zero nodes with velocity set by BC
	NodalVelBC *nextBC=firstVelocityBC;
	while(nextBC!=NULL)
		nextBC = nextBC->ZeroVelocityBC(mtime,passType);
	
	// on hold, exit
	if(holdAllVelocityBCs) return;
	
	// Now add all velocities to nodes with velocity BCs
	nextBC=firstVelocityBC;
	while(nextBC!=NULL)
		nextBC = nextBC->AddVelocityBC(mtime,passType);
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
