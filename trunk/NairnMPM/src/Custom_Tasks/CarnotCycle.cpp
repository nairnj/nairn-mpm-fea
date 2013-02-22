/********************************************************************************
	CarnotCycle.cpp
	NairnMPM

	Created by John Nairn on 2/11/13.
	Copyright (c) 2012 John A. Nairn, All rights reserved.
 
	Step 1. At thermal.reference, do isothermal expansion until
			average J = V1rel
    Step 2. Do adiabatic expansion until it cools to T2. Find V2rel
    Step 3. Do isothermal compression toV3rel. For ideal gas
              V3rel = V2rel/V1rel or specify V4rel
    Step 4. Do adiabatic compression until T=thermal.reference
 
	Vrel = 1 + exx for 1D compression in x direction
********************************************************************************/

#include "Custom_Tasks/CarnotCycle.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "Exceptions/CommonException.hpp"
#include "Exceptions/MPMTermination.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Materials/MaterialBase.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Global_Quantities/BodyForce.hpp"

#pragma mark Constructors and Destructors

// Constructors
CarnotCycle::CarnotCycle()
{
	V1rel = -1.;
	V2rel = 1.;
	V3rel = -1.;
	T2 = -1.;
}

// Return name of this task
const char *CarnotCycle::TaskName(void) { return "Carnot Cycle Simulation"; }

// Read task parameter - if pName is valid, set input for type
//    and return pointer to the class variable
char *CarnotCycle::InputParam(char *pName,int &input)
{
    if(strcmp(pName,"V1rel")==0)
    {	input=DOUBLE_NUM;
        return (char *)&V1rel;
    }
	
    else if(strcmp(pName,"V3rel")==0)
    {	input=DOUBLE_NUM;
        return (char *)&V3rel;
    }
	
    else if(strcmp(pName,"T2")==0)
    {	input=DOUBLE_NUM;
        return (char *)&T2;
    }
	
	// check remaining commands
    return CustomTask::InputParam(pName,input);
}

#pragma mark GENERIC TASK METHODS

// called once at start of MPM analysis - initialize and print info
CustomTask *CarnotCycle::Initialize(void)
{
    cout << "Carnot Cyle Simulation." << endl;
	
	// time interval
	cout << "   Step 1: Isothermal expansion to " << V1rel << "*V0" << endl;
	cout << "   Step 2: Adibatic expansion to temperature " << T2 << endl;
	if(V3rel>0.)
		cout << "   Step 3: Isothermal compression to " << V3rel << "*V0" << endl;
	else
		cout << "   Step 3: Isothermal compression to V3/V0 = V2/V1" << endl;
	cout << "   Step 2: Adibatic compression to temperature " << thermal.reference << endl;
	
	if(V1rel<1.)
		throw CommonException("The first expansion volume (V1rel) must be greater than 1","CarnotCycle::Initialize()");
	if(T2>thermal.reference || T2<0)
		throw CommonException("The second temperature, T2, must be less than the reference temperature.","CarnotCycle::Initialize()");
	if(V3rel>0. and V3rel<1.)
		throw CommonException("The thirds expansion volume (V3rel) must be greater than 1","CarnotCycle::Initialize()");
	
	// initial settings
	carnotStep = 0;
	ConductionTask::energyCoupling = FALSE;
	MaterialBase::isolatedSystemAndParticles = FALSE;
	
    return nextTask;
}

// Adjust time step now
CustomTask *CarnotCycle::StepCalculation(void)
{
    
    int p;
	MaterialBase *matID;
	
	double Vrel = 0.;
	int numgas = 0;
	double Tgas = 0.;
    
    // loop over material points
    for(p=0;p<nmpms;p++)
	{	// verify material is defined and sets if field number (in in multimaterial mode)
		matID=theMaterials[mpm[p]->MatID()];	// material object for this particle
		if(matID->Rigid()) continue;
		
        numgas++;
		Vrel += mpm[p]->GetRelativeVolume();
		Tgas += mpm[p]->pPreviousTemperature;
	}
	
	// get averages
	Vrel /= (double)numgas;
	Tgas /= (double)numgas;
	
	switch(carnotStep)
	{	case 0:
			carnotStep = 1;
			cout << "# Step 1: isothermal expansion" << endl;
            ConductionTask::energyCoupling = FALSE;
			break;
			
		case 1:
			// stop when Vrel reaches V1rel
			if(Vrel >= V1rel)
			{	// switch to adibatic expanion
				carnotStep = 2;
				ConductionTask::energyCoupling = TRUE;
				cout << "# Step 2: adibatic expansion" << endl;
			}
			break;
		
		case 2:
			// stop when T cools to T2
			if(Tgas<=T2)
			{	// switch to pause
				carnotStep = 3;
				ConductionTask::energyCoupling = FALSE;
				V2rel = Vrel;
				
				// find velocity
				for(p=0;p<nmpms;p++)
				{	// verify material is defined and sets if field number (in in multimaterial mode)
					matID=theMaterials[mpm[p]->MatID()];	// material object for this particle
					if(matID->Rigid())
						mpm[p]->ReverseParticle();
					else
						mpm[p]->StopParticle();
				}
				
				// estimate next volume
				if(V3rel<1.)
					V3rel = V2rel/V1rel;
				
				cout << "# Step 2: isothermal compression to " << V3rel << "*V0" << endl;
			}
			break;
		
		case 3:
			// stop when Vrel reaches V3rel
			if(Vrel <= V3rel)
			{	// switch to adibatic compression
				carnotStep = 4;
				ConductionTask::energyCoupling = TRUE;
				cout << "# Step 4: adibatic compression" << endl;
			}
			break;
		
		case 4:
			// done when T increases to thermal.reference
			if(Tgas>=thermal.reference)
				throw MPMTermination("Carnot cycle is complete","CarnotCycle::StepCalculation()");
			break;
		
		default:
			break;
	}

    return nextTask;
}

