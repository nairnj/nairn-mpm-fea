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
	xcentre= 0.0; 					// initial position of x centre EDIT
	ycentre = 0.0; 							// initial position of y centre EDIT
	zcentre= 12.25;					// depends for cs. COMMENT FOR HPT
	zDepth = 0.0;							// initial position of bottom of pin EDIT
	rpm=materialID=0.0;
	velSet.x=velSet.y=velSet.z=0.0;
	simTime=0.0;
	phi = 0.0;
}

// set the rpm from generators.cpp
void Rpm::SetRPM(double materialNumber, double revspermin)
{
	materialID=materialNumber;
	rpm=revspermin;
	//phi=(rpm/60)*2*PI_CONSTANT*timestep;		// for displacement rotation
	phi=(rpm/60)*2*PI_CONSTANT;					// for velocity rotation (negative for anticlockwise rotation)
}

bool Rpm::CheckRPM(int inputMatID)
{
	if((inputMatID)==materialID)
		return TRUE;
	else
		return FALSE;
}
/*
// not used as moved to SetRPM so only calculated once...
void Rpm::SetRotationAngle(double timestep)
{
	//Actually setting radial velocity now: 
	//This is multiplied by the radius of each point to give a velocity in x and y.
	//phi=(rpm/60)*2*PI_CONSTANT*timestep;		// for displacement rotation
	phi=(rpm/60)*2*PI_CONSTANT;					// for velocity rotation
	// show angle per step below 
	//double angleout = (phi*360)/(2*PI_CONSTANT);
	//cout << "*** Rotation angle per step = " << angleout <<" ***" << endl;
}
*/
/*
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
*/
void Rpm::AddRPM3(double &posz, double &posy, double &posx, Vector *vel, MPMBase *mptr, double timestep, double timer)
{
	sx=posx-xcentre;				// determine relative displacement
	sy=posy-ycentre;
	sz = posz - zcentre;
	
	//carry out initial rotation about y-axis, usually 2-3 degrees. And plunge, here only for away from plate.
 //START HERE FOR 3D ROTATION
	if(timer==0)
	{	//zcentre=12.5;
		
		//sz = posz - zcentre;
		double angle = 2*PI_CONSTANT*3/360;
		posx+=sx*cos(angle)-sz*sin(angle)-sx;					// apply rotation displacement to original position, negative phi for tilt in -ve x direction
		posz+=sx*sin(angle)+sz*cos(angle)-sz;
		
		//now apply the initial plunge if needed (outside of plate)
		posz-=0.04;
		
		ux = -sin(angle);
		uz = cos(angle);
	
	
	}
	
	// For continuous rotation about an axis tilted at 3 degrees:
	
	if(timer>0){ // start 3D rotation
	
	
	
	// can the centre point be on a different plane, NO.
	//Therefore each plane of particles must have it's own centre for each of x,y,z.
	//Can have one point accounting for the centre, and then a fixed displacement for each other plane.
	//Assume the centre is at the middle of the tool.
	//This will be difficult to differentiate between the different layers, so may use a few rigid bodies to model the tool
	
	
	//sx = 
	//sy = 
	//sz = 
	
	// This requires both a position and velocity update as divergence occurs in position with only velocity update.
	// rotates clockwise for positive phi. 
	posx += sx*(cos(phi*timestep)+ux*ux*(1-cos(phi*timestep))) + sy*(-uz*sin(phi*timestep)) + sz*(ux*uz*(1-cos(phi*timestep))) - sx;
	posy += sx*(uz*sin(phi*timestep)) + sy*cos(phi*timestep) + sz*(-ux*sin(phi*timestep)) - sy;
	posz += sx*(uz*ux*(1-cos(phi*timestep))) + sy*(ux*sin(phi*timestep)) + sz*(cos(phi*timestep)+uz*uz*(1-cos(phi*timestep))) - sz;
	
	// 3D EDIT
	//vel->x=0;						// zero velocities
	//vel->y=0;
	//vel->z=0;
	
	// positive phi means rotating clockwise. can change in SetRPM above or in input file (eg "rpm=-970")
	vel->x=-sy*uz*phi;					// velocity vector for rotation
	vel->y=sx*uz*phi - sz*ux*phi;
	//if (sy>=0)
		//vel->z-=sqrt(sy*sy + sx*sx) * phi;
	//else
		//vel->z+=sqrt(sy*sy + sx*sx) * phi;
	vel->z=sy*ux*phi;
	
	vel->x+=2.1167;
	posx+=timestep*2.1167;
	
	} //end 3D rotation
	//END HERE FOR 3D ROTATION
	
	
	
	//cout << "sx = " << sx << " sy = " << sy << endl;
	// untab the "EDIT" regions for 3D FSP...
	
	// 2D EDIT
	//vel->x=velSet.x;						// reset velocities
	//vel->y=velSet.y;


	// 3D EDIT
	//vel->x=0;						// zero velocities
	//vel->y=0;
	//vel->z=0;
	
	
	// EDIT FOR PLUNGE AND DWELL **************************
	//if(timer>timestep&&timer<=516*timestep) //create plunge stage 0.5 res
	//if(timer>timestep&&timer<=1031*timestep) //create plunge stage 1.0 res
	//if(timer>timestep&&timer<=103*timestep) //create plunge stage 1.0 res with ^5 density
	//if(timer>timestep&&timer<=103*timestep) // HPT plunge stage 1.0 res with ^5 density
	//if(timer>timestep&&timer<=5151*timestep) // FSP thermal testing
	//if(timer>timestep&&timer<=24.) //units are in mm/second copper exp1
		//vel->z=-0.1333; //HPT
	//	vel->z=-0.125; // FSP
			
									// difference in time between the two is dwell stage
	// EDIT FOR TRANSLATION *******************************
	//if(timer>=550*timestep)			// create translate stage
	//if(timer>=1100*timestep)			// create translate stage
	//if(timer>=35&&timer<90.)				// create translate stage
	//if(timer>=0)				// create translate stage
		//vel->x+=2.1167;
	
	//if(timer>=90.)				//create tool retraction stage
	//	vel->z=0.5;

/* // START here for HPT

	vel->x-=phi*sy;					// add new velocity vector for rotation
	vel->y+=phi*sx;
	
	
	//**************************************************
	// section to update particle position so to avoid lack of convergence:
	
	rotx=sx*cos(phi)-sy*sin(phi);	// new position relative to a zero origin
	roty=sx*sin(phi)+sy*cos(phi);
	
	// rotation in x-y
	posx+=sx*cos(phi*timestep)-sy*sin(phi*timestep)-sx;					// apply displacement to original position
	posy+=sx*sin(phi*timestep)+sy*cos(phi*timestep)-sy;
	
	//posx+=timestep*2.1167;												// TRANSLATION
	
	//cout << "vel x = " << vel->x << " vely = " << vel->y << endl;
*/ // END HERE FOR HPT

}



void Rpm::SetVel(Vector VelGen)
{
	velSet.x=VelGen.x;
	velSet.y=VelGen.y;
	velSet.z=VelGen.z;
}

void Rpm::UpdateCentre(double timestep)
{
	// 2D no plunge or dwell (tab out for 3D)
	//xcentre+=(timestep*velSet.x);
	//ycentre+=(timestep*velSet.y);
	
	//Keep in mind, currently with 3D sims, the rotating piece must start at (0,0)
	
	// 3D FSP with plunge or dwell (tab out for 2D)
	// must be calculated as the velocity in x is changed for translation stage. 
	// change times for translation stage above
	
	// EDIT FOR TRANSLATION AND HEAT FLUX *******************************
	//if(rotator->simTime>=550*timestep)
	//if(rotator->simTime>=1100*timestep)
	//if(rotator->simTime<=24.) // for HARD CODED HEAT FLUX EDIT
			//zDepth -=(timestep*0.125);
			
	//if(rotator->simTime>=35.&&rotator->simTime<90.) // for FSP TRANSLATION EDIT
	//if(rotator->simTime>=0) // for FSP TRANSLATION EDIT
 //START HERE FOR 3D ROTATION
		if(rotator->simTime>0)
			xcentre+=(timestep*2.1167);
		else   // if(rotator->simTime==0)	
		zcentre-=0.04;
	//END HERE FOR 3D ROTATION
	//cout << "Sim time = " << simTime << " x-centre = " << xcentre << " y-centre = " << ycentre << endl;
	//if(rotator->simTime>=90.) // for HARD CODED HEAT FLUX EDIT
		//	zDepth +=(timestep*0.5);

}

void Rpm::LogTime(double currentTime)
{
	simTime=currentTime;
}

// function for #hardcodedheat in ConductionTask.cpp
double Rpm::getxCentre(){return xcentre;}
double Rpm::getyCentre(){return ycentre;}
double Rpm::getzDepth(){return zDepth;}

