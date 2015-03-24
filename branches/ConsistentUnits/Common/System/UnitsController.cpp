/********************************************************************************
	UnitsController.cpp
	NairnFEA/NairnMPM

	Created by John Nairn on 3/15/15.
	Copyright (c) 2015 John A. Nairn, All rights reserved.
********************************************************************************/

#include "System/UnitsController.hpp"
#include "Read_XML/CommonReadHandler.hpp"

// class variables
int UnitsController::unitsType = LEGACY_UNITS;

// Constructor
UnitsController::UnitsController(void)
{
	
}

// When reading a variable, set scale factor (if in legacy units)
// and return the pointer
char *UnitsController::ScaledPtr(char *ptr,double &ptrScaling,double legacyScale)
{
	if(unitsType==LEGACY_UNITS) ptrScaling = legacyScale;
	return ptr;
}

// Return provided scalling factor or 1 if using consistent units
double UnitsController::Scaling(double legacyScale)
{	return unitsType==LEGACY_UNITS ? legacyScale : 1. ;
}

// Return label for a physical quantity
const char *UnitsController::Label(int type)
{
	switch(unitsType)
	{	case LEGACY_UNITS:
			switch(type)
			{	case ERR_UNITS:
					return "J/m^2";
				case CULENGTH_UNITS:
					return "mm";
				case CUMASS_UNITS:
					return "g";
				case TIME_UNITS:
					return "sec";
				case CUVELOCITY_UNITS:
					return "mm/sec";
				case TRACTIONSLOPE_UNITS:
					return "MPa/mm";
				case PRESSURE_UNITS:
					return "MPa";
				case HEATCAPACITY_UNITS:
					return "J/(kg-K)";
				case C0_UNITS:
					return "m/sec";
				case VISCOSITY_UNITS:
					return "cP";
				case CONDUCTIVITY_UNITS:
					return "W/(m-K)";
				case DENSITY_UNITS:
					return "g/cm^3";
				case STRESSINTENSITY_UNITS:
					return "MPa-sqrt(m)";
				case FEAFORCE_UNITS:
					// for FEA output
					return "N";
				case FEAWORK_UNITS:
					// for FEA output
					return "J";
				case INTERFACEPARAM_UNITS:
					return "MPa/mm";
				case TARGETKE_UNITS:
					// special case for feedback damping - other units same as work
					return "micro J";
				case BCARG_UNITS:
					return "ms/ms^-1";
				case BCTIME_UNITS:
					return "ms";
				case BCHEATFLUX_UNITS:
					return "W/m^2";
				case BCCONCFLUX_UNITS:
					return "kg/(m^2-s)";
				default:
					break;
			}
			break;
			
		default:
			break;
	}
	
	// no label gound
	return "????";
			
}

// When reading a file and find 'units' attribute, convert current
// units into those units (if possible)
double UnitsController::UnitsAttribute(char *value,int type)
{
	double attrScale = 1.;
	
	switch(type)
	{	case SEC_UNITS:
			// convert to seconds
			if(strcmp(value,"msec")==0)
				attrScale=1.e-3;
			else if(strcmp(value,"ms")==0)
				attrScale=1.e-3;
			else if(strcmp(value,"microsec")==0)
				attrScale=1.e-6;
			else if(strcmp(value,"us")==0)
				attrScale=1.e-6;
			break;
			
		case LENGTH_UNITS:
			// convert to mm
			if(strcmp(value,"m")==0)
				attrScale=1.e3;
			else if(strcmp(value,"cm")==0)
				attrScale=10.;
			else if(strcmp(value,"microns")==0)
				attrScale=1.e-3;
			else if(strcmp(value,"in")==0)
				attrScale=25.4;
			else if(strcmp(value,"ft")==0)
				attrScale=12.*25.4;
			break;
			
		case VELOCITY_UNITS:
			// convert to mm/sec
			if(strcmp(value,"m/sec")==0)
				attrScale=1.e3;
			else if(strcmp(value,"mm/msec")==0)
				attrScale=1.e3;
			else if(strcmp(value,"cm/sec")==0)
				attrScale=10.;
			else if(strcmp(value,"in/sec")==0)
				attrScale=25.4;
			else if(strcmp(value,"ft/sec")==0)
				attrScale=12.*25.4;
			break;
			
		case MASS_UNITS:
			// convert to g
			if(strcmp(value,"kg")==0)
				attrScale=1.e3;
			else if(strcmp(value,"mg")==0)
				attrScale=1.e-3;
			else if(strcmp(value,"lbs")==0)
				attrScale=453.594;
			else if(strcmp(value,"oz")==0)
				attrScale=28.3495;
			break;
			
		default:
			break;
	}
	
	return attrScale;
}

