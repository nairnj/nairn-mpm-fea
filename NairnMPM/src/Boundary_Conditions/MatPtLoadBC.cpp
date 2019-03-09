/********************************************************************************
    MatPtLoadBC.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Tues Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Boundary_Conditions/MatPtLoadBC.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Elements/ElementBase.hpp"
#include "System/UnitsController.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"
#include "Exceptions/CommonException.hpp"

// global
MatPtLoadBC *firstLoadedPt=NULL;

#pragma mark MatPtLoadBC::Constructors and Destructors

// Constructors
MatPtLoadBC::MatPtLoadBC(int num,int dof,int sty)
			: BoundaryCondition(sty,(double)0.,(double)0.)
{
    ptNum=num;
    direction=dof;		// note that 3 or 4 means z here, which may not be Z_DIRECTION (4) bit location
						// it is assumed z if not x (=1) or y (=2)
    holding=false;
}

// When particles reordered, fix point num if needed
// nonrigid particle p2+1 just moved to p1+1 and rigid particle
//		p1+1 moved to p2+1
MatPtLoadBC *MatPtLoadBC::ReorderPtNum(int p1, int p2)
{	if(ptNum==p2+1)
		ptNum = p1+1;
	else if(ptNum==p1+1)
		ptNum = p2+1;
	return (MatPtLoadBC *)GetNextObject();
}

#pragma mark MatPtLoadBC: Methods

// print to output
BoundaryCondition *MatPtLoadBC::PrintBC(ostream &os)
{
    char nline[200];
    
    sprintf(nline,"%7d %2d %2d %15.7e %15.7e",ptNum,direction,style,
			UnitsController::Scaling(1.e-6)*GetBCValueOut(),GetBCFirstTimeOut());
    os << nline;
	PrintFunction(os);
	
	// not allowed in rigid contact or BC particles
	MaterialBase *matref = theMaterials[mpm[ptNum-1]->MatID()];
	if(matref->IsRigid())
	{	if(!matref->IsRigidBlock())
		{	throw CommonException("Cannot set external force on rigid contact or rigid boundary condition particles",
								  "MatPtLoadBC::PrintBC");
		}
		else if(style==SILENT)
		{	throw CommonException("Cannot set silent external force on rigid block particles",
								  "MatPtLoadBC::PrintBC");
		}
	}
	
    // Initial value is F or F/numParticles if net, but if function
    //      initial value is s or s/numParticles if net
	// (in Legacy, F in uN and s=1e6, otherwise, F unscaled and s=1)
	
	// rescale ... for fuction only using s or s/numParticles if net force
	if(style==FUNCTION_VALUE)
		scale = GetBCValue();
	
    return (BoundaryCondition *)GetNextObject();
}

// set external load to zero in both directions
MatPtLoadBC *MatPtLoadBC::ZeroMPLoad(void)
{
	ZeroVector(mpm[ptNum-1]->GetPFext());
    return (MatPtLoadBC *)GetNextObject();
}

// increment external load on a particle
// input is analysis time in seconds
MatPtLoadBC *MatPtLoadBC::AddMPLoad(double bctime)
{
	if(style!=SILENT)
	{	Vector *pFext=mpm[ptNum-1]->GetPFext();
		
		if(direction==X_DIRECTION)
			 pFext->x+=BCValue(bctime);
		else if(direction==Y_DIRECTION)
			 pFext->y+=BCValue(bctime);
		else
			 pFext->z+=BCValue(bctime);
	}
	else
	{	MPMBase *mptr = mpm[ptNum-1];
		double mp=mptr->mp;												// in g
		int matnum=mptr->MatID();
		double cd=theMaterials[matnum]->CurrentWaveSpeed(FALSE,mpm[ptNum-1],0);	// in mm/sec (2D or isotropic only)
		double cs=theMaterials[matnum]->ShearWaveSpeed(FALSE,mpm[ptNum-1],0);		// in mm/sec
		Vector pVel=mptr->vel;
		Vector *pFext=mptr->GetPFext();
		
		// get forces in g-mm/sec^2
		if(direction==X_DIRECTION)
		{	pFext->x-=mp*cd*pVel.x/(2.*mptr->GetParticleXSize());
			pFext->y-=mp*cs*pVel.y/(2.*mptr->GetParticleXSize());
			pFext->z-=mp*cs*pVel.z/(2.*mptr->GetParticleXSize());
		}
		else if(direction==Y_DIRECTION)
		{	pFext->x-=mp*cs*pVel.x/(2.*mptr->GetParticleYSize());
			pFext->y-=mp*cd*pVel.y/(2.*mptr->GetParticleYSize());
			pFext->z-=mp*cs*pVel.z/(2.*mptr->GetParticleYSize());
		}
		else
		{	pFext->x-=mp*cs*pVel.x/(2.*mptr->GetParticleZSize());
			pFext->y-=mp*cs*pVel.y/(2.*mptr->GetParticleZSize());
			pFext->z-=mp*cd*pVel.z/(2.*mptr->GetParticleZSize());
		}
	}
		
    return (MatPtLoadBC *)GetNextObject();
}

// reverse active linear loads. Leave rest alone
// input is analysis time in seconds
// if LINEAR_VALUE, set finalTime to time when load returns to zero
MatPtLoadBC *MatPtLoadBC::ReverseLinearLoad(double bctime,double *finalTime,bool holdFirst)
{
    switch(style)
    {	case LINEAR_VALUE:
            if(bctime>=GetBCFirstTime())
			{   SetBCOffset(BCValue(bctime));
                if(holdFirst)
				{	// change to offset (set above) with zero slope to get a constant load
                    holdValue = GetBCValue();
                    holding = true;
                    SetBCValue(0.);
                }
                else
				{	// get unloading slope (units/time)
					double bcvalue = holding ? -holdValue : -GetBCValue() ;
					
					// change to offset+(slope)*(time-ftime) by setting offset (above), new slope, and new first time 
					SetBCValueCU(bcvalue);
                    SetBCFirstTimeCU(bctime);
                    
					// find time when BC returns to zero
                    // new BC is offset+value*(mstime-ftime), which is zero when mstime = ftime-offset/value
                    *finalTime = GetBCFirstTime() - GetBCOffset()/GetBCValue();
                }
                
            }
            break;
        default:
            break;
    }
    return (MatPtLoadBC *)GetNextObject();
}

// hold active linear loads at current values. Leave rest alone
// input is analysis time in seconds
MatPtLoadBC *MatPtLoadBC::MakeConstantLoad(double bctime)
{
    switch(style)
    {	case LINEAR_VALUE:
            if(bctime>=GetBCFirstTime())
            {	style=CONSTANT_VALUE;
				SetBCValueCU(BCValue(bctime));
            }
            break;
        default:
            break;
    }
    return (MatPtLoadBC *)GetNextObject();
}

#pragma mark MatPtLoadBC:Flux method in sub classes

// add "flux" condition - must override in sub class flux BC
MatPtLoadBC *MatPtLoadBC::AddMPFluxBC(double bctime)
{	return (MatPtLoadBC *)GetNextObject();
}


#pragma mark MatPtLoadBC:Accessors

// get current position of particle
void MatPtLoadBC::GetPosition(unordered_map<string,double> &vars)
{	vars["x"] = mpm[ptNum-1]->pos.x;
	vars["y"] = mpm[ptNum-1]->pos.y;
	vars["z"] = mpm[ptNum-1]->pos.z;
	vars["q"] = mpm[ptNum-1]->GetParticleRotationZ();
}

// set value (and scale legacy N to uN, and MPa to Pa)
void MatPtLoadBC::SetBCValue(double bcvalue)
{	BoundaryCondition::SetBCValue(UnitsController::Scaling(1.e6)*bcvalue);
}

#pragma mark MatPtLoadBC: Class Methods

// Calculate forces applied to particles at time stepTime in g-mm/sec^2
void MatPtLoadBC::SetParticleFext(double stepTime)
{
	// zero all with external force BC
    MatPtLoadBC *nextLoad=firstLoadedPt;
    while(nextLoad!=NULL)
    	nextLoad=nextLoad->ZeroMPLoad();
	
	// add external force BCs
    nextLoad=firstLoadedPt;
    while(nextLoad!=NULL)
    	nextLoad=nextLoad->AddMPLoad(stepTime);
}

