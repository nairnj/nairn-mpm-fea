/********************************************************************************
    ThermalRamp.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Jan 08 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.    
********************************************************************************/

#include "stdafx.h"
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
	reference=0.;
}

/*******************************************************************
	Methods
*******************************************************************/

// print thermal ramp info to output file
void ThermalRamp::Output(void)
{
	char hline[200];
	
	// always print reference temperature
	sprintf(hline,"Reference Temperature: %g K",reference);
	cout << hline << endl;
}





