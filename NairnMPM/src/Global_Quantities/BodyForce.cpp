/********************************************************************************
    BodyForce.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Jan 08 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.    
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Global_Quantities/BodyForce.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Nodes/NodalPoint.hpp"
#include "System/UnitsController.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Read_XML/Expression.hpp"

extern double timestep;

// Single global object
// Handles body forces (now only gravity and damping)
BodyForce bodyFrc;

#pragma mark BodyForce:Initialize

// Constructors
BodyForce::BodyForce()
{
	gravity = false;
	hasGridBodyForce = false;
	
	damping=0.;                 // constant damping
	useDamping = false;
	dampingCoefficient=0.;		// 1/Q in Nose-Hoover feedback
	useFeedback = false;
	alpha=0.;					// evolving damping coefficient
    maxAlpha=-1.;               // max evolving alpha
	function=NULL;              // target kinetic energy function
	gridfunction = NULL;        // for constant damping that depends on time
    
	pdamping=0.;                // same for particle damping
	usePDamping = false;
	pdampingCoefficient=0.;
	usePFeedback = false;
	
    palpha=0.;                  // evoloving particle damping coefficient
    maxPAlpha=-1.;              // max evolving particle alpha
	pfunction=NULL;              // target kinetic energy function
	pgridfunction=NULL;
    
	useGridFeedback=TRUE;		// base feedback on grid kinetic energy
								// provide option to change to allow particle kintic energy instead
	
	XPICOrder = 0;				// default to FLIP
	isUsingVstar = 0;			// XPIC being used (always 0 in NairnMPM)
	xpicVectors = 1;			// to store vk and one is added to store pk
	usingFMPM = false;
	gridBCOption = GRIDBC_COMBINED;
	
	gridBodyForceFunction[0]=NULL;
	gridBodyForceFunction[1]=NULL;
	gridBodyForceFunction[2]=NULL;
}

// Destructor (and it is virtual)
BodyForce::~BodyForce()
{	if(function!=NULL) delete function;
    if(pfunction!=NULL) delete pfunction;
	if(gridfunction!=NULL) delete gridfunction;
    if(pgridfunction!=NULL) delete pgridfunction;
	if(gridBodyForceFunction[0]!=NULL) delete gridBodyForceFunction[0];
	if(gridBodyForceFunction[1]!=NULL) delete gridBodyForceFunction[1];
	if(gridBodyForceFunction[2]!=NULL) delete gridBodyForceFunction[2];
}

#pragma mark BodyForce:Gravity and Body Force Functions

// turn gravity on (initially zero forces)
void BodyForce::Activate(void)
{
    // zero to start, components get set directly by BodyXForce, etc., XML commands
    ZeroVector(&gforce);
	gravity = true;
}

// Get gravity (constant body force) and grid based functions (in mm/sec^2)
// Multiplied by mass when added to the material velocity field
// If call in parallel, use try block to catch possible function errors
void BodyForce::GetGridBodyForce(Vector *theFrc,Vector *fpos,double utime)
{
	// constant body force
	if(gravity)
	{	theFrc->x = gforce.x;
		theFrc->y = gforce.y;
		theFrc->z = gforce.z;
	}
	else
		theFrc->x = theFrc->y = theFrc->z = 0.;
	
	// exit if no grid functions
	if(!hasGridBodyForce) return;

    // set variables (see Expression vmap)
	double vars[5];
	vars[0] = 4.5;
	vars[1] = utime*UnitsController::Scaling(1.e3);		//t
	vars[2] = fpos->x;		//x
	vars[3] = fpos->y;		//y
	vars[4] = fpos->z;		//z
	
	// body force functions
	if(gridBodyForceFunction[0]!=NULL)
		theFrc->x += gridBodyForceFunction[0]->EvaluateFunction(vars);
	if(gridBodyForceFunction[1]!=NULL)
		theFrc->y += gridBodyForceFunction[1]->EvaluateFunction(vars);
	if(gridBodyForceFunction[2]!=NULL)
		theFrc->z += gridBodyForceFunction[2]->EvaluateFunction(vars);
}


// set function for grid body force
// throws std::bad_alloc, SAXException()
void BodyForce::SetGridBodyForceFunction(char *bcFunction,int input)
{
	if(bcFunction==NULL)
	{	ThrowSAXException("Grid body force function is missing");
		return;
	}
	if(strlen(bcFunction)==0)
	{	ThrowSAXException("Grid body force function is missing");
		return;
	}
	
	// turn it on
	hasGridBodyForce = true;			// true if any functions provided
	Expression *newFunction =  Expression::CreateExpression(bcFunction,"Grid body force function is not valid");
	
	// assign to direction
	if(input==GRID_X_BODY_FORCE_FUNCTION_BLOCK)
	{	if(gridBodyForceFunction[0]!=NULL)
        	ThrowSAXException("Duplicate grid body x force function");
		gridBodyForceFunction[0] = newFunction;
	}
	else if(input==GRID_Y_BODY_FORCE_FUNCTION_BLOCK)
	{	if(gridBodyForceFunction[1]!=NULL)
        	ThrowSAXException("Duplicate grid body y force function");
		gridBodyForceFunction[1] = newFunction;
	}
	else if(input==GRID_Z_BODY_FORCE_FUNCTION_BLOCK)
	{	if(gridBodyForceFunction[2]!=NULL)
        	ThrowSAXException("Duplicate grid body z force function");
		gridBodyForceFunction[2] = newFunction;
	}
}

#pragma mark BodyForce:Grid and Particle Damping

// Return true if need to add any forces in the post forces task
bool BodyForce::HasGridDampingForces()
{
	// True if body forces (gravity or functions)
	if(gravity || bodyFrc.hasGridBodyForce) return true;
	
	// nothing needed to change forces
	return false;
}

// Get sum of damping and feedback damping to apply to the grid
double BodyForce::GetGridDamping(double utime)
{
	if(!useDamping) return 0.;
	
	double totalDamping;
	
	// simple damping
	if(gridfunction==NULL)
		totalDamping = damping;
	else
	{	totalDamping =  gridfunction->TValue(utime*UnitsController::Scaling(1000.));
	}
	
	// add feedback damping
	if(useFeedback) totalDamping += alpha;
	
	return totalDamping;
}

// Get sum of damping and feedback damping to apply in particle velocity update
double BodyForce::GetParticleDamping(double utime)
{
	if(!usePDamping) return 0.;
	
	double totalDamping;
	
	// simple damping
	if(pgridfunction==NULL)
		totalDamping = pdamping;
	else
	{	totalDamping =  pgridfunction->TValue(utime*UnitsController::Scaling(1000.));
	}
	
	// add feedback damping
	if(usePFeedback) totalDamping += palpha;
	
	return totalDamping;
}

// display gravity settings
void BodyForce::Output(void)
{
	char hline[200];
	size_t hsize=200;
	
    // Gravity
	if(gravity)
	{	snprintf(hline,hsize,"Body force per %s: (%g,%g,%g) %s/%s^2",UnitsController::Label(CUMASS_UNITS),
				gforce.x,gforce.y,gforce.z,UnitsController::Label(CULENGTH_UNITS),UnitsController::Label(TIME_UNITS));
		cout << hline << endl;
	}

	// Body force functions
	if(gridBodyForceFunction[0]!=NULL)
	{	const char *expr = gridBodyForceFunction[0]->GetString();
		cout << "Grid body x force per " << UnitsController::Label(CUMASS_UNITS) << ": " << expr
					<< " " << UnitsController::Label(CULENGTH_UNITS) << "/" << UnitsController::Label(TIME_UNITS) << "^2" << endl;
	}
	if(gridBodyForceFunction[1]!=NULL)
	{	const char *expr = gridBodyForceFunction[1]->GetString();
		cout << "Grid body y force per " << UnitsController::Label(CUMASS_UNITS) << ": " << expr
					<< " " << UnitsController::Label(CULENGTH_UNITS) << "/" << UnitsController::Label(TIME_UNITS) << "^2" << endl;
	}
	if(gridBodyForceFunction[2]!=NULL)
	{	const char *expr = gridBodyForceFunction[2]->GetString();
		cout << "Grid body z force per " << UnitsController::Label(CUMASS_UNITS) << ": " << expr
					<< " " << UnitsController::Label(CULENGTH_UNITS) << "/" << UnitsController::Label(TIME_UNITS) << "^2" << endl;
	}
	
	// Grid damping
	if(gridfunction!=NULL)
	{	const char *expr = gridfunction->GetString();
		cout << "Grid damping = " << expr << " /" << UnitsController::Label(TIME_UNITS) << endl;
	}
	else if(damping!=0.)
	{	snprintf(hline,hsize,"Grid damping: %g /%s",damping,UnitsController::Label(TIME_UNITS));
		cout << hline << endl;
	}
    else if(!useFeedback)
    {   // turn off if nothing above and no feedback too
        useDamping = false;
    }

    // Grid feedback damping
	if(useFeedback)
	{	snprintf(hline,hsize,"Grid feedback damping with coefficient: %g /%s^2",dampingCoefficient,UnitsController::Label(CULENGTH_UNITS));
		cout << hline << endl;
		if(function!=NULL)
		{	const char *expr = function->GetString();
			cout << "   Target kinetic energy = " << expr << " " << UnitsController::Label(TARGETKE_UNITS) << endl;
		}
        else
            cout << "   Target kinetic energy = 0" << endl;
        if(maxAlpha>0.)
        {	snprintf(hline,hsize,"   Maximum grid damping alpha: %g /%s",maxAlpha,UnitsController::Label(TIME_UNITS));
            cout << hline << endl;
        }
	}
    
    // Particle damping
	if(pgridfunction!=NULL)
	{	const char *expr = pgridfunction->GetString();
		cout << "Particle damping = " << expr << " /" << UnitsController::Label(TIME_UNITS) << endl;
	}
	else if(pdamping!=0.)
	{	snprintf(hline,hsize,"Particle damping: %g /%s",pdamping,UnitsController::Label(TIME_UNITS));
		cout << hline << endl;
	}
    else if(!usePFeedback)
    {   // turn off if nothing above and no feedback too
        usePDamping = false;
    }
    
    // Particle feedback damping
	if(usePFeedback)
	{	snprintf(hline,hsize,"Particle feedback damping with coefficient: %g /%s^2",pdampingCoefficient,UnitsController::Label(CULENGTH_UNITS));
		cout << hline << endl;
		if(pfunction!=NULL)
		{	const char *expr = pfunction->GetString();
			cout << "   Target kinetic energy = " << expr << " " << UnitsController::Label(TARGETKE_UNITS) << endl;
		}
        else
            cout << "   Target kinetic energy = 0" << endl;
        if(maxPAlpha>0.)
        {	snprintf(hline,hsize,"   Maximum particle damping alpha: %g /%s",maxPAlpha,UnitsController::Label(TIME_UNITS));
            cout << hline << endl;
        }
	}
}

// update alpha normalized to number of particles
// only call when feedback is known to be on
void BodyForce::UpdateAlpha(double delTime,double utime)
{
    // if neither use feedback, then exit
	if(!useFeedback && !usePFeedback) return;
	
	// get total kinetic energy (same for both types)
	double kineticEnergy=0.;
	double totalMass=0.;
	if(useGridFeedback)
	{	int i;
		for(i=1;i<=nnodes;i++)
			nd[i]->AddKineticEnergyAndMass(kineticEnergy,totalMass);
	}
	else
	{	int p;
		for(p=1;p<nmpmsNR;p++)
		{	if(mpm[p]->InReservoir()) continue;		// don't add those in resevoir
			kineticEnergy += mpm[p]->KineticEnergy();
			totalMass += mpm[p]->mp;
		}
	}
	
    // Grid damping
    if(useFeedback)
    {   // target kinetic energy (Legacy use time in ms and targetEnergy converted to nJ)
        double targetEnergy;
        if(function!=NULL)
        {	targetEnergy =  function->TValue(utime*UnitsController::Scaling(1000.))*UnitsController::Scaling(1000.);
        }
        else
            targetEnergy = 0.;

        // compare target to kinetic energy
        // factor 2 is for historic reasons to match old damping factor
        alpha += 2.*dampingCoefficient*(kineticEnergy-targetEnergy)*delTime/totalMass;
        if(alpha<0.)
            alpha=0.;
        else if(maxAlpha>0)
        {   if(alpha>maxAlpha)
                alpha=maxAlpha;
        }
    }

    // Particle damping
    if(usePFeedback)
    {   // target kinetic energy (Legacy use time in ms and targetEnergy converted to nJ)
        double targetEnergy;
        if(pfunction!=NULL)
        {	targetEnergy =  pfunction->TValue(utime*UnitsController::Scaling(1000.))*UnitsController::Scaling(1000.);
        }
        else
            targetEnergy=0.;
        
        // compare target to kinetic energy
        // factor 2 is for historic reasons to match old damping factor
        palpha += 2.*pdampingCoefficient*(kineticEnergy-targetEnergy)*delTime/totalMass;
        if(palpha<0.)
            palpha=0.;
        else if(maxPAlpha>0)
        {   if(palpha>maxPAlpha)
                palpha=maxPAlpha;
        }
    }
}

// set target function for feedback damping
// throws std::bad_alloc, SAXException()
void BodyForce::SetTargetFunction(char *bcFunction,bool gridDamp)
{   
	if(bcFunction==NULL)
	{	ThrowSAXException("Target energy function of time is missing");
		return;
	}
	if(strlen(bcFunction)==0)
	{	ThrowSAXException("Target energy function of time is missing");
		return;
	}
	if((gridDamp && function!=NULL) || (!gridDamp && pfunction!=NULL))
	{	ThrowSAXException("Duplicate target energy function of time for grid or particle damping");
		return;
	}

	// create the function and check it
	if(gridDamp)
	{	function =  Expression::CreateExpression(bcFunction,
							"Target energy function of time for grid damping is not valid");
	}
	else
	{	pfunction =  Expression::CreateExpression(bcFunction,
							"Target energy function of time for particle damping is not valid");
	}
}

// set maximum alpha (units 1/sec)
// throws SAXException()
void BodyForce::SetMaxAlpha(double theMax,bool gridDamp)
{   if(theMax<=0.)
	{	ThrowSAXException("Maximum feedback damping alpha must be positive");
		return;
	}
    if(gridDamp)
        maxAlpha = theMax;
    else
        maxPAlpha = theMax;
}

// set function for grid damping
// throws std::bad_alloc, SAXException()
void BodyForce::SetGridDampingFunction(char *bcFunction,bool gridDamp)
{
	if(bcFunction==NULL)
	{	ThrowSAXException("Grid or particle damping function of time is missing");
		return;
	}
	if(strlen(bcFunction)==0)
	{	ThrowSAXException("Grid or particle damping function of time is missing");
		return;
	}
	if((gridDamp && gridfunction!=NULL) || (!gridDamp && pgridfunction!=NULL))
	{	ThrowSAXException("Duplicate grid or particle damping function of time");
		return;
	}

	// create the function and check it
	if(gridDamp)
	{	gridfunction =  Expression::CreateExpression(bcFunction,
									"Grid damping function of time is not valid");
	}
	else
	{	pgridfunction =  Expression::CreateExpression(bcFunction,
									"Particle damping function of time is not valid");
	}
}

// Implement extended PIC (but only allow 1 or higher)
void BodyForce::SetXPICOrder(int newOrder) { XPICOrder = newOrder>=0 ? newOrder : 0 ; }

// return XPIC order (1 is normal PIC, 2 or higher is XPIC/FMPM )
int BodyForce::GetXPICOrder(void) { return XPICOrder; }

// return if simulations is using high-order XPIC/FMPM (2 or higher)
// called during set up and to allocate XPIC tasks
int BodyForce::UsingVstar(void) { return isUsingVstar; }
void BodyForce::SetUsingVstar(int setting) { isUsingVstar = setting; }

// change the number needed if add PeriodicXPIC custom task
int BodyForce::XPICVectors(void) { return xpicVectors; }
void BodyForce::SetXPICVectors(int vnum) { xpicVectors = vnum; }

// Is this XPIC or FMPM (note that not changed for intersperse FLIP steps)
bool BodyForce::UsingFMPM(void) { return usingFMPM; }
void BodyForce::SetUsingFMPM(bool useit) { usingFMPM = useit; }

// FMPM(k>1) and XPIC(k>1) grid BC option
int BodyForce::GridBCOption(void) { return gridBCOption; }
void BodyForce::SetGridBCOption(int newOption) { gridBCOption = newOption; }
