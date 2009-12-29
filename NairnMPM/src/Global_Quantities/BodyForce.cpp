/********************************************************************************
    BodyForce.cpp
    NairnMPM
    
    Created by John Nairn on Thu Jan 08 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.    
********************************************************************************/

#include "Global_Quantities/BodyForce.hpp"
#include "Read_XML/CommonReadHandler.hpp"
#include "Read_XML/mathexpr.hpp"

// Single global object
// Handles body forces (now only gravity and damping)
BodyForce bodyFrc;

// global expression variables
double BodyForce::varTime=0.;
PRVar keTimeArray[1] = { NULL };

/*******************************************************************
	ThermalRamp: Constructors and Destructors
*******************************************************************/

// Constructors
BodyForce::BodyForce()
{
	gravity=FALSE;
	damping=0.;
	useFeedback=FALSE;
	dampingCoefficient=0.;		// 1/Q in Nose-Hoover feedback
	alpha=0.;					// evolving damping coefficient
}

/*******************************************************************
	Methods
*******************************************************************/

// turn gravity on (initially zero forces)
void BodyForce::Activate(void)
{
	gforcex=0.;
	gforcey=0.;
	gravity=TRUE;
	function=NULL;
}

// Destructor (and it is virtual)
BodyForce::~BodyForce()
{	if(function!=NULL) delete function;
}

// If gravity return TRUE and current forces
short BodyForce::GetGravity(double *xfrc,double *yfrc)
{
	if(!gravity) return FALSE;
	*xfrc=gforcex;
	*yfrc=gforcey;
	return TRUE;
}

// the damping
double BodyForce::GetDamping(void) { return damping; }

// display gravity settings
void BodyForce::Output(void)
{
	char hline[200];
	
	if(gravity)
	{	sprintf(hline,"Body force per g: (%e,%e) mm/sec^2",gforcex,gforcey);
		cout << hline << endl;
	}
	
	sprintf(hline,"Grid damping: %g /sec",damping);
    cout << hline << endl;
	
	if(useFeedback)
	{	sprintf(hline,"Feedback damping with coefficient: %g /mm^2",dampingCoefficient);
		cout << hline << endl;
		if(function!=NULL)
		{	char *expr=function->Expr('#');
			cout << "   Target kinetic energy = " << expr << " micro J" << endl;
			delete [] expr;
		}
	}

}

// return the damping coefficient
double BodyForce::GetAlpha(void) { return alpha; }

// initialize to zero
void BodyForce::TrackAlpha(void)
{	kineticEnergy=0.;
	totalMass=0.;
}

// increment if has damping using grid velocity extrapolated to the particle
void BodyForce::TrackAlpha(Vector *vstar,double mp)
{
	if(!useFeedback) return;
	
	// twice particle kinetic energy in g mm^2/sec^2
	kineticEnergy+=mp*(vstar->x*vstar->x+vstar->y*vstar->y);
	totalMass+=mp;
}

// update alpha normalized to number of particles
void BodyForce::UpdateAlpha(double delTime,double utime)
{
	if(!useFeedback) return;
	
	// target kinetic energy in micro J
	double targetEnergy;
	if(function!=NULL)
	{	varTime=1000.*utime;
		targetEnergy=function->Val();
	}
	else
		targetEnergy=0.;

	// actual kinetic energy in micro J is kineticEnergy*0.5e-3
	// this damping factor has units of 1/mm^2 and extra factor of 2.e3 to make same
	//    magnitude as previous damping method
	alpha+=dampingCoefficient*(kineticEnergy-2.e3*targetEnergy)*delTime/totalMass;
	if(alpha<0.) alpha=0.;
}

// set target function for feedback damping
void BodyForce::SetTargetFunction(char *bcFunction)
{
	if(bcFunction==NULL)
		throw SAXException("Target energy function of time is missing");
	if(strlen(bcFunction)==0)
		throw SAXException("Target energy function of time is missing");
	if(function!=NULL)
		throw SAXException("Duplicate target energy functions of time");
	
	// create variable
	if(keTimeArray[0]==NULL)
	{	keTimeArray[0]=new RVar("t",&varTime);
	}
		
	// create the function and check it
	function=new ROperation(bcFunction,1,keTimeArray);
	if(function->HasError())
		throw SAXException("Target energy function of time is not valid");
}


