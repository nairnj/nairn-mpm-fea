/********************************************************************************
    ThermalRamp.cpp
    nairn-mpm-fea
 
	Used to have thermal ramp, which is now done by CustomThermalRamp
	This class only stores reference temperature in a global variable
    
    Created by John Nairn on Thu Jan 08 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.    
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Global_Quantities/ThermalRamp.hpp"
#include "System/UnitsController.hpp"

#include "MPM_Classes/MPMBase.hpp"
#include <iostream>
#include <cstddef>

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
    size_t hsize=200;
	
	// always prints reference temperature
	snprintf(hline,hsize,"Reference Temperature: %g K",reference);
	cout << hline << endl;
}





