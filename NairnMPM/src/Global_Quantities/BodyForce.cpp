/********************************************************************************
    BodyForce.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Jan 08 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.    
********************************************************************************/

#include "stdafx.h"
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
	gravity=FALSE;
	hasGridBodyForce=FALSE;
	
	damping=0.;                 // constant damping
	useDamping = false;
	dampingCoefficient=0.;		// 1/Q in Nose-Hoover feedback
	useFeedback=FALSE;
	alpha=0.;					// evolving damping coefficient
    maxAlpha=-1.;               // max evolving alpha
	function=NULL;              // target kinetic energy function
	gridfunction = NULL;        // for constant damping that depends on time
    
	pdamping=0.;                // same for particle damping
	usePDamping = false;
	pdampingCoefficient=0.;
	usePFeedback=FALSE;
    palpha=0.;                  // evoloving particle damping coefficient
    maxPAlpha=-1.;              // max evolving particle alpha
	pfunction=NULL;              // target kinetic energy function
	pgridfunction=NULL;
    
	useGridFeedback=TRUE;		// base feedback on grid kinetic energy
								// provide option to change to allow particle kintic energy instead
	
	fractionPIC = 0.;			// fraction PIC implemeted by damping
	usePICDamping = false;		// off by default
	XPICOrder = 1;				// default to normal PIC (if PIC is used) (always 1 in NairnMPM)
	
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
	gravity=TRUE;
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

	unordered_map<string,double> vars;
	vars["t"] = utime*UnitsController::Scaling(1000.);
	vars["x"] = fpos->x;
	vars["y"] = fpos->y;
	vars["z"] = fpos->z;
	
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
	hasGridBodyForce = TRUE;
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

// Get sum of damping and feedback damping to apply to the grid
// (Don't call from parallel code due to function)
double BodyForce::GetDamping(double utime)
{
    double totalDamping;
    
    // simple damping
	if(gridfunction==NULL)
        totalDamping = damping;
    else
	{	totalDamping =  gridfunction->TValue(utime*UnitsController::Scaling(1000.));
    }
    
    // add feedback damping
    if(useFeedback) totalDamping += alpha;
	
    // PIC Damping
	if(usePICDamping)
	{	totalDamping -= (double)XPICOrder*fractionPIC/timestep;
	}
    
    return totalDamping;
}

// Get sum of damping and feedback damping to apply in particle velocity update
double BodyForce::GetParticleDamping(double utime)
{
    double totalDamping;
    
    // simple damping
	if(pgridfunction==NULL)
        totalDamping = pdamping;
    else
	{	totalDamping =  pgridfunction->TValue(utime*UnitsController::Scaling(1000.));
    }
    
    // add feedback damping
    if(usePFeedback) totalDamping += palpha;
    
    // PIC Damping
	if(usePICDamping) totalDamping += fractionPIC/timestep;
	
    return totalDamping;
}

// Get grid damping without the PIC term
// (Don't call from parallel code due to function in GetDamping())
double BodyForce::GetNonPICDamping(double utime)
{
    bool hold = usePICDamping;
    usePICDamping = false;
    double nonPICDamping = GetDamping(utime);
    usePICDamping = hold;
    return nonPICDamping;
}

// Get particle damping without the PIC term
double BodyForce::GetNonPICParticleDamping(double utime)
{
    bool hold = usePICDamping;
    usePICDamping = false;
    double nonPICDamping = GetParticleDamping(utime);
    usePICDamping = hold;
    return nonPICDamping;
}

// Get PIC damping term alone
double BodyForce::GetPICDamping(void)
{   return usePICDamping ? fractionPIC/timestep : 0.;
}

// display gravity settings
void BodyForce::Output(void)
{
	char hline[200];
	
    // Gravity
	if(gravity)
	{	sprintf(hline,"Body force per %s: (%g,%g,%g) %s/%s^2",UnitsController::Label(CUMASS_UNITS),
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
	{	sprintf(hline,"Grid damping: %g /%s",damping,UnitsController::Label(TIME_UNITS));
		cout << hline << endl;
	}
    else if(!useFeedback)
    {   // turn off if nothing above and no feedback too
        useDamping = false;
    }
    
    // Grid feedback damping
	if(useFeedback)
	{	sprintf(hline,"Grid feedback damping with coefficient: %g /%s^2",dampingCoefficient,UnitsController::Label(CULENGTH_UNITS));
		cout << hline << endl;
		if(function!=NULL)
		{	const char *expr = function->GetString();
			cout << "   Target kinetic energy = " << expr << " " << UnitsController::Label(TARGETKE_UNITS) << endl;
		}
        else
            cout << "   Target kinetic energy = 0" << endl;
        if(maxAlpha>0.)
        {	sprintf(hline,"   Maximum grid damping alpha: %g /%s",maxAlpha,UnitsController::Label(TIME_UNITS));
            cout << hline << endl;
        }
	}
    
    // Particle damping
	if(pgridfunction!=NULL)
	{	const char *expr = pgridfunction->GetString();
		cout << "Particle damping = " << expr << " /" << UnitsController::Label(TIME_UNITS) << endl;
	}
	else if(pdamping!=0.)
	{	sprintf(hline,"Particle damping: %g /%s",pdamping,UnitsController::Label(TIME_UNITS));
		cout << hline << endl;
	}
    else if(!usePFeedback)
    {   // turn off if nothing above and no feedback too
        usePDamping = false;
    }
    
    // Particle feedback damping
	if(usePFeedback)
	{	sprintf(hline,"Particle feedback damping with coefficient: %g /%s^2",pdampingCoefficient,UnitsController::Label(CULENGTH_UNITS));
		cout << hline << endl;
		if(pfunction!=NULL)
		{	const char *expr = pfunction->GetString();
			cout << "   Target kinetic energy = " << expr << " " << UnitsController::Label(TARGETKE_UNITS) << endl;
		}
        else
            cout << "   Target kinetic energy = 0" << endl;
        if(maxPAlpha>0.)
        {	sprintf(hline,"   Maximum particle damping alpha: %g /%s",maxPAlpha,UnitsController::Label(TIME_UNITS));
            cout << hline << endl;
        }
	}

    // PIC damping
	if(usePICDamping)
	{	cout << "PIC damping fraction: " << fractionPIC;
		cout << endl;
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
		{	kineticEnergy += mpm[p]->KineticEnergy();
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

// Implement fraction PIC by setting damping values
// throws SAXException()
void BodyForce::SetFractionPIC(double fract)
{
	if(fract<0. || fract>1.)
	{	ThrowSAXException("Fraction PIC damping parameter must be from 0 to 1.");
		return;
	}
	
	if(fract>0.)
	{	fractionPIC = fract;
		usePICDamping = true;
	}
	else
	{	fractionPIC = 0.;
		usePICDamping = false;
	}
}

// bool if fracturePIC is ot zero
bool BodyForce::IsUsingPICDamping(void) { return usePICDamping;}

// get fraction PIC
double BodyForce::GetFractionPIC(void) { return fractionPIC; }

