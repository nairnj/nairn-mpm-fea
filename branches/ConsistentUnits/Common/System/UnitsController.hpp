/********************************************************************************
	UnitsController.hpp
	NairnFEA/NairnMPM

	Created by John Nairn on 3/15/15.
	Copyright (c) 2015 John A. Nairn, All rights reserved.
********************************************************************************/

#ifndef _UNITSCONTROLLER_
#define _UNITSCONTROLLER_

enum { LEGACY_UNITS=0 };

enum { ERR_UNITS=0,OUTLENGTH_UNITS,TRACTIONSLOPE_UNITS,PRESSURE_UNITS,HEATCAPACITY_UNITS,
		C0_UNITS,VISCOSITY_UNITS,CONDUCTIVITY_UNITS,DENSITY_UNITS,STRESSINTENSITY_UNITS,
		FORCE_UNITS, WORK_UNITS, INTERFACEPARAM_UNITS };

class UnitsController
{
	public:
	
		//  Constructors and Destructor
		UnitsController();
	
		// class methods
		static char *ScaledPtr(char *,double &,double);
		static double Scaling(double);
		static const char *Label(int);

	
	private:
		static int unitsType;
		
};

#endif
