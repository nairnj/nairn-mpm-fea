// Rigid Body rotation

/********************************************************************************
    RPM.cpp
    NairnMPM
    
    Created by Tim Fagan on Thu Dec 1 2011.
********************************************************************************/

#include "Read_MPM/RPM.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Nodes/NodalPoint.hpp"

Rpm *rotator;
bool Rpm::rpmApplied=FALSE;


//Constructors and Destructors

Rpm::Rpm() 
{
	rpmApplied=TRUE;
	simTime=0.0;
	xcentre=ycentre=zcentre=0.0;
	rpm=materialID=0.0;
	velSet.x=velSet.y=velSet.z=0.0;
	simTime=0.0;
}

// set the rpm from generators.cpp
void Rpm::SetRPM(double materialNumber, double revspermin)
{
	materialID=materialNumber;
	rpm=revspermin;
	//phi=(rpm/60)*2*PI_CONSTANT*timestep;		// for displacement rotation
	phi=(rpm/60)*2*PI_CONSTANT;					// for velocity rotation
}

bool Rpm::CheckRPM(int inputMatID)
{
	if((inputMatID)==materialID)
		return TRUE;
	else
		return FALSE;
}

// not used as moved to SetRPM so only calculated once...
void Rpm::SetRotationAngle(double timestep)
{
	/*Actually setting radial velocity now: 
	This is multiplied by the radius of each point to give a velocity in x and y.*/
	//phi=(rpm/60)*2*PI_CONSTANT*timestep;		// for displacement rotation
	phi=(rpm/60)*2*PI_CONSTANT;					// for velocity rotation
	// show angle per step below 
	//double angleout = (phi*360)/(2*PI_CONSTANT);
	//cout << "*** Rotation angle per step = " << angleout <<" ***" << endl;
}


void Rpm::AddRPM(double &posy, double &posx, double &extwork, Vector pFext)
{
	//phi=(rpm/60)*2*PI_CONSTANT*(timestep);
	sx=posx-xcentre;				// determine relative position to a zero origin
	sy=posy-ycentre;
	rotx=sx*cos(phi)-sy*sin(phi);	// new position relative to a zero origin
	roty=sx*sin(phi)+sy*cos(phi);
	posx+=rotx-sx;					// apply displacement to original position
	posy+=roty-sy;
	extwork+=rotx*pFext.x+roty*pFext.y;

}

void Rpm::AddRPM2(double &posy, double &posx, double &extwork, Vector &pFext, Vector *vel, MPMBase *mptr, double timestep, double timer)
{
	sx=posx-xcentre;				// determine relative displacement
	sy=posy-ycentre;

	// untab the "EDIT" regions for 3D FSP...
	
	// 2D
	vel->x=velSet.x;						// reset velocities
	vel->y=velSet.y;

	
	// 3D
	//vel->x=0;						// zero velocities
	//vel->y=0;
	//vel->z=0;
	
	
	// EDIT FOR PLUNGE AND DWELL **************************
	//if(timer>timestep&&timer<=4*timestep) //create plunge stage
	//		vel->z=-400e5;
									// difference in time between the two is dwell stage
	// EDIT FOR TRANSLATION *******************************
	//if(timer>=7*timestep)			// create translate stage
	//		vel->x=20e5;
	
	vel->x-=phi*sy;					// add new velocity vector for rotation
	vel->y+=phi*sx;

}

void Rpm::AddRPM3(double &posy, double &posx, Vector *vel, MPMBase *mptr, double timestep, double timer)
{
	sx=posx-xcentre;				// determine relative displacement
	sy=posy-ycentre;

	// untab the "EDIT" regions for 3D FSP...
	
	// 2D
	vel->x=velSet.x;						// reset velocities
	vel->y=velSet.y;

	
	// 3D
	//vel->x=0;						// zero velocities
	//vel->y=0;
	//vel->z=0;
	
	
	// EDIT FOR PLUNGE AND DWELL **************************
	//if(timer>timestep&&timer<=4*timestep) //create plunge stage
	//		vel->z=-400e5;
									// difference in time between the two is dwell stage
	// EDIT FOR TRANSLATION *******************************
	//if(timer>=7*timestep)			// create translate stage
	//		vel->x=20e5;
	
	vel->x-=phi*sy;					// add new velocity vector for rotation
	vel->y+=phi*sx;

}



void Rpm::SetVel(Vector VelGen)
{
	velSet.x=VelGen.x;
	velSet.y=VelGen.y;
	velSet.z=VelGen.z;
}

void Rpm::UpdateCentre(double timestep)
{
	// 2D no plunge or dwell
	xcentre+=(timestep*velSet.x);
	ycentre+=(timestep*velSet.y);
	
	// 3D FSP with plunge or dwell
	// must be calculated as the velocity in x is changed for translation stage. 
	// change times for translation stage above
	// EDIT FOR TRANSLATION *******************************
	//if(rotator->simTime>=7*timestep)
	//		xcentre+=(timestep*20e5);
}

void Rpm::LogTime(double currentTime)
{
	simTime=currentTime;
}



