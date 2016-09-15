/********************************************************************************
    ThermalRamp.hpp
    nairn-mpm-fea
    
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
		double reference;		// reference or stress free temperature
	
        // constructors and destructors
        ThermalRamp();
    
        // methods
		void Output(void);
	
};

extern ThermalRamp thermal;

#endif
