/********************************************************************************
    BodyForce.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Jan 08 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.    
********************************************************************************/

#include "Global_Quantities/BodyForce.hpp"
#include "Read_XML/mathexpr.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Nodes/NodalPoint.hpp"

extern double timestep;

// Single global object
// Handles body forces (now only gravity and damping)
BodyForce bodyFrc;

// global expression variables
double BodyForce::varTime=0.;
PRVar keTimeArray[1] = { NULL };
PRVar gTimeArray[1] = { NULL };

#pragma mark BodyForce:Initialize

// Constructors
BodyForce::BodyForce()
{
	gravity=FALSE;
	
	damping=0.;                 // constant damping
	useDamping=FALSE;
	dampingCoefficient=0.;		// 1/Q in Nose-Hoover feedback
	useFeedback=FALSE;
	alpha=0.;					// evolving damping coefficient
    maxAlpha=-1.;               // max evolving alpha
	function=NULL;              // target kinetic energy function
	gridfunction=NULL;          // for constant damping that depends on time
    
	pdamping=0.;                // same for particle damping
	usePDamping=FALSE;
	pdampingCoefficient=0.;
	usePFeedback=FALSE;
    palpha=0.;                  // evoloving particle damping coefficient
    maxPAlpha=-1.;              // max evolving particle alpha
	pgridfunction=NULL;
    
	useGridFeedback=TRUE;		// base feedback on grid kinetic energy
                                // provide option to change to allow particle kintic energy instead
	
	fractionPIC = 0.;			// fraction PIC implemeted by damping
	usePICDamping = false;		// off by default
    
}

// Destructor (and it is virtual)
BodyForce::~BodyForce()
{	if(function!=NULL) delete function;
    if(pfunction!=NULL) delete pfunction;
	if(gridfunction!=NULL) delete gridfunction;
    if(pgridfunction!=NULL) delete pgridfunction;
}

#pragma mark BodyForce:Gravity and Body Force Functions

// turn gravity on (initially zero forces)
void BodyForce::Activate(void)
{
    // zero to start, components get set directly by BodyXForce, etc., XML commands
    ZeroVector(&gforce);
	gravity=TRUE;
}

// If gravity add to material point force buffer
void BodyForce::AddGravity(Vector *theFrc,double mp,double wt)
{
	if(!gravity) return;
	
	double gscale = mp*wt;
	theFrc->x += gscale*gforce.x;
	theFrc->y += gscale*gforce.y;
	theFrc->z += gscale*gforce.z;
}

#pragma mark BodyForce:Grid and Particle Damping

// Get sum of damping and feedback damping to apply to the grid
double BodyForce::GetDamping(double utime)
{
    double totalDamping;
    
    // simple damping
	if(gridfunction==NULL)
        totalDamping = damping;
    else
	{   varTime=1000.*utime;
        totalDamping = gridfunction->Val();
    }
    
    // add feedback damping
    if(useFeedback) totalDamping += alpha;
	
    // PIC damping
	if(usePICDamping) totalDamping -= fractionPIC/timestep;
    
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
	{   varTime=1000.*utime;
        totalDamping = pgridfunction->Val();
    }
    
    // add feedback damping
    if(usePFeedback) totalDamping += palpha;
    
    // PIC Damping
	if(usePICDamping) totalDamping += fractionPIC/timestep;
	
    return totalDamping;
}

// Get grid damping without the PIC term
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

// display gravity settings
void BodyForce::Output(void)
{
	char hline[200];
	
    // Gravity
	if(gravity)
	{	sprintf(hline,"Body force per g: (%g,%g,%g) mm/sec^2",gforce.x,gforce.y,gforce.z);
		cout << hline << endl;
	}
	
    // Grid damping
	if(gridfunction!=NULL)
	{	char *expr=gridfunction->Expr('#');
		cout << "Grid damping = " << expr << " /sec" << endl;
		delete [] expr;
	}
	else
	{	sprintf(hline,"Grid damping: %g /sec",damping);
		cout << hline << endl;
	}
	if(useFeedback)
	{	sprintf(hline,"Grid feedback damping with coefficient: %g /mm^2",dampingCoefficient);
		cout << hline << endl;
		if(function!=NULL)
		{	char *expr=function->Expr('#');
			cout << "   Target kinetic energy = " << expr << " micro J" << endl;
			delete [] expr;
		}
        else
            cout << "   Target kinetic energy = 0" << endl;
        if(maxAlpha>0.)
        {	sprintf(hline,"   Maximum grid damping alpha: %g /sec",maxAlpha);
            cout << hline << endl;
        }
	}
    
    // Particle damping
	if(pgridfunction!=NULL)
	{	char *expr=pgridfunction->Expr('#');
		cout << "Particle damping = " << expr << " /sec" << endl;
		delete [] expr;
	}
	else
	{	sprintf(hline,"Particle damping: %g /sec",pdamping);
		cout << hline << endl;
	}
	if(usePFeedback)
	{	sprintf(hline,"Particle feedback damping with coefficient: %g /mm^2",pdampingCoefficient);
		cout << hline << endl;
		if(pfunction!=NULL)
		{	char *expr=pfunction->Expr('#');
			cout << "   Target kinetic energy = " << expr << " micro J" << endl;
			delete [] expr;
		}
        else
            cout << "   Target kinetic energy = 0" << endl;
        if(maxPAlpha>0.)
        {	sprintf(hline,"   Maximum particle damping alpha: %g /sec",maxPAlpha);
            cout << hline << endl;
        }
	}

	if(usePICDamping)
	{	cout << "PIC damping fraction: " << fractionPIC << endl;
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
    {   // target kinetic energy in micro J
        double targetEnergy;
        if(function!=NULL)
        {	varTime=1000.*utime;
            targetEnergy=function->Val();
        }
        else
            targetEnergy=0.;

        // actual kinetic energy in micro J is kineticEnergy*1.0e-3
        // for target energy in g mm^2/sec^2 is 1000*targetEnergy
        // this damping factor has units of 1/mm^2 and extra factor of 2.e3 to make same
        //    magnitude as previous damping method
        alpha += 2.*dampingCoefficient*(kineticEnergy-1000.*targetEnergy)*delTime/totalMass;
        if(alpha<0.)
            alpha=0.;
        else if(maxAlpha>0)
        {   if(alpha>maxAlpha)
                alpha=maxAlpha;
        }
    }

    // Particle damping
    if(usePFeedback)
    {   // target kinetic energy in micro J
        double targetEnergy;
        if(pfunction!=NULL)
        {	varTime=1000.*utime;
            targetEnergy=pfunction->Val();
        }
        else
            targetEnergy=0.;
        
        // actual kinetic energy in micro J is kineticEnergy*1.0e-3
        // for target energy in g mm^2/sec^2 is 1000*targetEnergy
        // this damping factor has units of 1/mm^2 and extra factor of 2.e3 to make same
        //    magnitude as previous damping method
        palpha += 2.*pdampingCoefficient*(kineticEnergy-1000.*targetEnergy)*delTime/totalMass;
        if(palpha<0.)
            palpha=0.;
        else if(maxPAlpha>0)
        {   if(palpha>maxPAlpha)
                palpha=maxPAlpha;
        }
    }
}

// set target function for feedback damping
void BodyForce::SetTargetFunction(char *bcFunction,bool gridDamp)
{   
	if(bcFunction==NULL)
		ThrowSAXException("Target energy function of time is missing");
	if(strlen(bcFunction)==0)
		ThrowSAXException("Target energy function of time is missing");
	if((gridDamp && function!=NULL) || (!gridDamp && pfunction!=NULL))
		ThrowSAXException("Duplicate target energy function of time for grid or particle damping");
	
	// create variable
	if(keTimeArray[0]==NULL)
	{	keTimeArray[0]=new RVar("t",&varTime);
	}

	// create the function and check it
    if(gridDamp)
    {   function=new ROperation(bcFunction,1,keTimeArray);
        if(function->HasError())
            ThrowSAXException("Target energy function of time for grid damping is not valid");
    }
    else
    {   pfunction=new ROperation(bcFunction,1,keTimeArray);
        if(pfunction->HasError())
            ThrowSAXException("Target energy function of time for particle is not valid");
    }
}

// set maximum alpha (units 1/sec)
void BodyForce::SetMaxAlpha(double theMax,bool gridDamp)
{   if(theMax<=0.)
        ThrowSAXException("Maximum feedback damping alpha must be positive");
    if(gridDamp)
        maxAlpha = theMax;
    else
        maxPAlpha = theMax;
}

// set function for grid damping
void BodyForce::SetGridDampingFunction(char *bcFunction,bool gridDamp)
{
	if(bcFunction==NULL)
		ThrowSAXException("Grid or particle damping function of time is missing");
	if(strlen(bcFunction)==0)
		ThrowSAXException("Grid or particle damping function of time is missing");
	if((gridDamp && gridfunction!=NULL) || (!gridDamp && pgridfunction!=NULL))
		ThrowSAXException("Duplicate grid or particle damping function of time");
	
	// create variable
	if(gTimeArray[0]==NULL)
	{	gTimeArray[0]=new RVar("t",&varTime);
	}
	
	// create the function and check it
    if(gridDamp)
	{   gridfunction=new ROperation(bcFunction,1,gTimeArray);
        if(gridfunction->HasError())
            ThrowSAXException("Grid damping function of time is not valid");
    }
    else
	{   pgridfunction=new ROperation(bcFunction,1,gTimeArray);
        if(pgridfunction->HasError())
            ThrowSAXException("Particle damping function of time is not valid");
    }
}

// Implement fraction PIC by setting damping values
void BodyForce::SetFractionPIC(double fract)
{
	if(fract<0. || fract>1.)
		ThrowSAXException("Fraction PIC damping parameter must be from 0 to 1.");
	
	if(fract>0.)
	{	fractionPIC = fract;
		usePICDamping = true;
	}
	else
		usePICDamping = false;
}


