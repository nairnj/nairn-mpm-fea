/********************************************************************************
    ThermalRamp.hpp
    NairnMPM
    
    Created by John Nairn on Thu Jan 08 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.  
	
	Header for ramping up constant temperature differences
		Global thermal

	Dependencies
		none
********************************************************************************/

#ifndef _THERMALRAMP_

#define _THERMALRAMP_

class ThermalRamp
{
    public:
		// constants (not changed in MPM time step, but set in XML commands)
		double isoDeltaT;		// final temperature change
		double isoRampTime;		// time to reach temperature change
		double reference;		// reference or stress free temperature
		double rampStart;		// when to start the ramp
		
        // constructors and destructors
        ThermalRamp();
    
        // methods
		void UpdateParticleTemperature(double *,double);
		void Activate(void);
        bool Active(void);
		void SetParameters(double);
		void CheckDone(double);
		void Output(void);
	
	private:
		// variables (changed in MPM time step)
		short isoRamp;			// ramps up thermal load
		double stepDTemp;		// total temperature change in this step
	
		// constants (not changed in MPM time step)
		double isoTempRate;		// rate when there is a ramp
};

extern ThermalRamp thermal;

#endif
