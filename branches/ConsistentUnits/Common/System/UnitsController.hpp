/********************************************************************************
	UnitsController.hpp
	NairnFEA/NairnMPM

	Created by John Nairn on 3/15/15.
	Copyright (c) 2015 John A. Nairn, All rights reserved.
********************************************************************************/

#ifndef _UNITSCONTROLLER_
#define _UNITSCONTROLLER_

enum { LEGACY_UNITS=0 };

enum { ERR_UNITS=0,LENGTH_UNITS,TRACTIONSLOPE_UNITS,PRESSURE_UNITS };

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
