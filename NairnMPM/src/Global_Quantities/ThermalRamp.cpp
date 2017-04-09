/********************************************************************************
    ThermalRamp.cpp
    nairn-mpm-fea
 
	Used to have thermal ramp, which is no done by CustomThermalRamp
	This class only stores reference temperature in a global variable
    
    Created by John Nairn on Thu Jan 08 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.    
********************************************************************************/

#include "stdafx.h"
#include "Global_Quantities/ThermalRamp.hpp"
#include "System/UnitsController.hpp"

// Single thermal ramp global object
ThermalRamp thermal;

#pragma mark ThermalRamp: Constructors and Destructors

// Constructors
ThermalRamp::ThermalRamp()
{
	reference=0.;
}

#pragma mark ThermalRamp: Methods

// print thermal ramp info to output file
void ThermalRamp::Output(void)
{
	char hline[200];
	
	// always prints reference temperature
	sprintf(hline,"Reference Temperature: %g K",reference);
	cout << hline << endl;
}





