/********************************************************************************
    MatPtLoadBC.cpp
    NairnMPM
    
    Created by John Nairn on Tues Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "Boundary_Conditions/MatPtLoadBC.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Materials/MaterialBase.hpp"
#include "NairnMPM_Class/MeshInfo.hpp"

// global
MatPtLoadBC *firstLoadedPt=NULL;

#pragma mark MatPtLoadBC::Constructors and Destructors

// Constructors
MatPtLoadBC::MatPtLoadBC(long num,int dof,int sty)
			: BoundaryCondition(sty,(double)0.,(double)0.)
{
    ptNum=num;
    direction=dof;		// note that 3 or 4 means z here, which may not be Z_DIRECTION (4) bit location
}

#pragma mark MatPtLoadBC: Methods

// print to output
BoundaryCondition *MatPtLoadBC::PrintBC(ostream &os)
{
    char nline[200];
    
    sprintf(nline,"%5ld %2d %2d %15.7e %15.7e",ptNum,direction,style,value,ftime);
    os << nline;
	PrintFunction(os);
	
	// rescale
	value*=1.e6;		// Multiply by 1e6 to get N (kg-m/sec^2) to g-mm/sec^2
	if(style==FUNCTION_VALUE)
		scale=value;		// ... value is 1.e6/numParticles 
	else
		scale*=1.e6;		// ... same in case using a function
	
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
	{	double mstime=1000.*bctime;
		Vector *pFext=mpm[ptNum-1]->GetPFext();
		
		if(direction==X_DIRECTION)
			 pFext->x+=BCValue(mstime);
		else if(direction==Y_DIRECTION)
			 pFext->y+=BCValue(mstime);
		else
			 pFext->z+=BCValue(mstime);
	}
	else
	{	double mp=mpm[ptNum-1]->mp;								// in g
		int matnum=mpm[ptNum-1]->MatID();
		double cd=theMaterials[matnum]->WaveSpeed();			// in mm/sec
		double cs=theMaterials[matnum]->ShearWaveSpeed();		// in mm/sec
		Vector pVel=mpm[ptNum-1]->vel;
		Vector *pFext=mpm[ptNum-1]->GetPFext();
		
		// get forces in g-mm/sec^2
		if(direction==X_DIRECTION)
		{	pFext->x-=mp*cd*pVel.x/mpmgrid.partx;
			pFext->y-=mp*cs*pVel.y/mpmgrid.partx;
			pFext->z-=mp*cs*pVel.z/mpmgrid.partx;
		}
		else if(direction==Y_DIRECTION)
		{	pFext->x-=mp*cs*pVel.x/mpmgrid.party;
			pFext->y-=mp*cd*pVel.y/mpmgrid.party;
			pFext->z-=mp*cs*pVel.z/mpmgrid.party;
		}
		else
		{	pFext->x-=mp*cs*pVel.x/mpmgrid.partz;
			pFext->y-=mp*cs*pVel.y/mpmgrid.partz;
			pFext->z-=mp*cd*pVel.z/mpmgrid.partz;
		}
	}
		
    return (MatPtLoadBC *)GetNextObject();
}

// reverse active linear loads. Leave rest alone
// input is analysis time in seconds
MatPtLoadBC *MatPtLoadBC::ReverseLinearLoad(double bctime)
{
	double mstime=1000.*bctime;
	
    switch(style)
    {	case LINEAR_VALUE:
            if(mstime>=ftime)
			{   offset=BCValue(mstime);
                value=-value;
                ftime=mstime;
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
	double mstime=1000.*bctime;
	
    switch(style)
    {	case LINEAR_VALUE:
            if(mstime>=ftime)
            {	style=CONSTANT_VALUE;
                value=BCValue(mstime);
            }
            break;
        default:
            break;
    }
    return (MatPtLoadBC *)GetNextObject();
}

#pragma mark MatPtLoadBC: Class Methods

// Calculate forces applied to particles at time stepTime in g-mm/sec^2
void MatPtLoadBC::SetParticleFext(double stepTime)
{
    MatPtLoadBC *nextLoad=firstLoadedPt;
    while(nextLoad!=NULL)
    	nextLoad=nextLoad->ZeroMPLoad();
    nextLoad=firstLoadedPt;
    while(nextLoad!=NULL)
    	nextLoad=nextLoad->AddMPLoad(stepTime);
}

#pragma mark MatPtLoadBC: Accessors

// get current position of particle
void MatPtLoadBC::GetPosition(double *xpos,double *ypos,double *zpos,double *rot)
{	*xpos=mpm[ptNum-1]->pos.x;
	*ypos=mpm[ptNum-1]->pos.y;
	*zpos=mpm[ptNum-1]->pos.z;
	*rot=mpm[ptNum-1]->GetRotation();
}


