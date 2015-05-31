/********************************************************************************
    ThermalRamp.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Jan 08 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.    
********************************************************************************/

#include "Global_Quantities/ThermalRamp.hpp"
#include "System/UnitsController.hpp"

// Single thermal ramp global object
ThermalRamp thermal;

/*******************************************************************
	ThermalRamp: Constructors and Destructors
*******************************************************************/

// Constructors
ThermalRamp::ThermalRamp()
{
	isoRamp=FALSE;
	reference=0.;
	isoDeltaT=0.;
	isoRampTime=-1.;
	rampStart=0.;
	stepDTemp=0.;
}

/*******************************************************************
	Methods
*******************************************************************/

// update particle temperature if ramp is on
void ThermalRamp::UpdateParticleTemperature(double *pTemperature,double delTime)
{	if(isoRamp) *pTemperature+=isoTempRate*delTime;
}

// turn ramp on (unless already on)
void ThermalRamp::Activate(void) { isoRamp=TRUE; }
bool ThermalRamp::Active(void) { return isoRamp; }

// set parameters for ramp at start of MPM analysis
void ThermalRamp::SetParameters(double delTime)
{
    if(isoRamp)
	{	int nsteps,nstart;
	
		// set rate and time when ramp is done
    	if(isoRampTime>0.)
        {   // apply over isoRampTime starting at rampTime
			nsteps=(int)(isoRampTime/delTime);
			if(nsteps<1) nsteps=1;
        }
        else
        {   // apply in one time step
			nsteps=1;
        }
		
		// set rate, start time, and ramp time
		isoTempRate=isoDeltaT/((double)nsteps*delTime);
		nstart=(int)(rampStart/delTime);
		rampStart=((double)nstart-0.1)*delTime;
		isoRampTime=rampStart+(double)nsteps*delTime;
    }
}

// turn off ramp when done and zero total increment for this step
void ThermalRamp::CheckDone(double currentTime)
{
    if(isoRamp)
    {	if(currentTime>isoRampTime)
			isoRamp=FALSE;
		stepDTemp=0.;
    }
}

// print thermal ramp info to output file
void ThermalRamp::Output(void)
{
	char hline[200];
	
	// always print reference temperature
	sprintf(hline,"Reference Temperature: %g K",reference);
	cout << hline << endl;
	
	if(!isoRamp) return;
	
	sprintf(hline,"Initial isothermal temperature difference: %g C",isoDeltaT);
	cout << hline << endl;
	if(isoRampTime>0.)
	{   sprintf(hline,"Ramped between %g and %g %s",rampStart*UnitsController::Scaling(1.e3),
				(rampStart+isoRampTime)*UnitsController::Scaling(1.e3),UnitsController::Label(ALTTIME_UNITS));
		cout << hline << endl;
	}
	else if(rampStart>0.)
	{   sprintf(hline,"Applied in one step at %g %s",rampStart*UnitsController::Scaling(1.e3),UnitsController::Label(ALTTIME_UNITS));
		cout << hline << endl;
	}
	else
		cout << "Applied in first time step" << endl;
}





