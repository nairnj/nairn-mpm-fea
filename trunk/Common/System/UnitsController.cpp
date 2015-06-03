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
int UnitsController::lengthExp=0,UnitsController::massExp=0,UnitsController::timeExp=0;
char UnitsController::length[3],UnitsController::mass[3],UnitsController::timeu[3];
char UnitsController::velocity[6],UnitsController::force[10],UnitsController::pressure[10];
char UnitsController::energy[10],UnitsController::density[10];
char UnitsController::errRate[15],UnitsController::tractionSlope[15],UnitsController::viscosity[15];
char UnitsController::heatCapacity[18],UnitsController::power[15],UnitsController::conductivity[20],UnitsController::convection[20];
char UnitsController::heatFlux[15],UnitsController::stressIntensity[20],UnitsController::bcarg[10],UnitsController::concFlux[20];

// Constructor
UnitsController::UnitsController(void)
{
	
}

// output units at start of the file
void UnitsController::OutputUnits(void)
{
	if(unitsType==LEGACY_UNITS)
		cout << "Legacy";
	else
		cout << length << "-" << mass << "-" << timeu;
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
			{	case CULENGTH_UNITS:
					return "mm";
				case CUMASS_UNITS:
					return "g";
				case TIME_UNITS:
					return "sec";
				case ALTTIME_UNITS:
					return "ms";
				case DENSITY_UNITS:
					return "g/cm^3";
				case CUVELOCITY_UNITS:
					return "mm/sec";
				case C0_UNITS:
					return "m/sec";
				case ALTVELOCITY_UNITS:
					return "m/sec";
				case FEAFORCE_UNITS:
					// for FEA output
					return "N";
				case PRESSURE_UNITS:
					return "MPa";
				case FEAWORK_UNITS:
					// for FEA output
					return "J";
				case TARGETKE_UNITS:
					// special case for feedback damping - other units same as work
					return "micro J";
				case ERR_UNITS:
					return "J/m^2";
				case STRESSINTENSITY_UNITS:
					return "MPa-sqrt(m)";
				case TRACTIONSLOPE_UNITS:
					return "MPa/mm";
				case INTERFACEPARAM_UNITS:
					return "MPa/mm";
				case VISCOSITY_UNITS:
					return "cP";
				case HEATCAPACITY_UNITS:
					return "J/(kg-K)";
				case CONDUCTIVITY_UNITS:
					return "W/(m-K)";
				case CONVECTION_UNITS:
					return "W/(m^2-K)";
				case BCHEATFLUX_UNITS:
					return "W/m^2";
				case BCARG_UNITS:
					return "ms/ms^-1";
				case BCCONCFLUX_UNITS:
					return "kg/(m^2-s)";
				default:
					break;
			}
			break;
			
		default:
			switch(type)
			{	case CULENGTH_UNITS:
					return length;
				case CUMASS_UNITS:
					return mass;
				case TIME_UNITS:
				case ALTTIME_UNITS:
					return timeu;
				case DENSITY_UNITS:
					return density;
				case CUVELOCITY_UNITS:
				case ALTVELOCITY_UNITS:
				case C0_UNITS:
					return velocity;
				case FEAFORCE_UNITS:
					return force;
				case PRESSURE_UNITS:
					return pressure;
				case FEAWORK_UNITS:
				case TARGETKE_UNITS:
					return energy;
				case ERR_UNITS:
					return errRate;
				case STRESSINTENSITY_UNITS:
					return stressIntensity;
				case TRACTIONSLOPE_UNITS:
				case INTERFACEPARAM_UNITS:
					return tractionSlope;
				case VISCOSITY_UNITS:
					return viscosity;
				case HEATCAPACITY_UNITS:
					return heatCapacity;
				case CONDUCTIVITY_UNITS:
					return conductivity;
				case CONVECTION_UNITS:
					return convection;
				case BCHEATFLUX_UNITS:
					return heatFlux;
				case BCARG_UNITS:
					return bcarg;
				case BCCONCFLUX_UNITS:
					return concFlux;
					
				default:
					break;
			}
			break;
	}
	
	// no label gound
	return "????";
			
}

// Create consistent units
bool UnitsController::SetConsistentUnits(char *len,char *ms,char *tm)
{
	if(unitsType != LEGACY_UNITS) return false;
	
	unitsType = CONSISTENT_UNITS;
	
	// length
	if(strcmp(len,"km")==0)
		lengthExp = 3;
	else if(strcmp(len,"m")==0)
		lengthExp = 0;
	else if(strcmp(len,"dm")==0)
		lengthExp = -1;
	else if(strcmp(len,"cm")==0)
		lengthExp = -2;
	else if(strcmp(len,"mm")==0)
		lengthExp = -3;
	else if(strcmp(len,"um")==0)
		lengthExp = -6;
	else if(strcmp(len,"microns")==0)
	{	strcpy(len,"um");
		lengthExp = -6;
	}
	else if(strcmp(len,"nm")==0)
		lengthExp = -9;
	else if(strcmp(len,"L")==0)
	{	lengthExp = 0;
		unitsType = USER_UNITS;
	}
	else
		return false;
	strcpy(length,len);
	
	// mass
	if(strcmp(ms,"kg")==0)
		massExp = 0;
	else if(strcmp(ms,"g")==0)
		massExp = -3;
	else if(strcmp(ms,"mg")==0)
		massExp = -6;
	else if(strcmp(ms,"ug")==0)
		massExp = -9;
	else if(strcmp(ms,"M")==0)
	{	massExp = 0;
		unitsType = USER_UNITS;
	}
	else
		return false;
	strcpy(mass,ms);
	
	// time
	if(strcmp(tm,"s")==0)
		timeExp = 0;
	else if(strcmp(tm,"sec")==0)
	{	tm[1]=0;
		timeExp = 0;
	}
	else if(strcmp(tm,"ms")==0)
		timeExp = -3;
	else if(strcmp(tm,"msec")==0)
	{	tm[2]=0;
		timeExp = -3;
	}
	else if(strcmp(tm,"us")==0)
		timeExp = -6;
	else if(strcmp(ms,"T")==0)
	{	timeExp = 0;
		unitsType = USER_UNITS;
	}
	else
		return false;
	strcpy(timeu,tm);
	
	// build labels
	
	// density
	strcpy(density,mass);
	strcat(density,"/");
	strcat(density,length);
	strcat(density,"^3");
	
	// velocity
	strcpy(velocity,len);
	strcat(velocity,"/");
	strcat(velocity,timeu);
	
	// force
	int fexp = massExp+lengthExp-2*timeExp;
	if(unitsType==USER_UNITS)
		strcpy(force,"F");
	else
	{	GetPrefix(fexp,"N",force);
		if(strlen(force)==0)
		{	GetPrefix(fexp+5,"dyne",force);
			if(strlen(force)==0)
				sprintf(force,"10^%dN",fexp);
		}
	}
	
	// pressure
	int pexp = fexp-2*lengthExp;
	if(unitsType==USER_UNITS)
		strcpy(pressure,"F/L^2");
	else
	{	GetPrefix(pexp,"Pa",pressure);
		if(strlen(pressure)==0)
		{	GetPrefix(pexp+1,"Ba",pressure);
			if(strlen(pressure)==0)
				sprintf(pressure,"10^%dPa",pexp);
		}
	}
	
	// energy
	int enerexp = fexp+lengthExp;
	if(unitsType==USER_UNITS)
		strcpy(energy,"F-L");
	else
	{	GetPrefix(enerexp,"J",energy);
		if(strlen(energy)==0)
		{	GetPrefix(enerexp+7,"erg",energy);
			if(strlen(energy)==0)
				sprintf(energy,"10^%dJ",enerexp);
		}
	}
	
	// energy release rate
	strcpy(errRate,energy);
	strcat(errRate,"/");
	strcat(errRate,length);
	strcat(errRate,"^2");
	
	// stress intensity
	strcpy(stressIntensity,pressure);
	strcat(stressIntensity,"-sqrt(");
	strcat(stressIntensity,length);
	strcat(stressIntensity,")");

	// traction slope
	strcpy(tractionSlope,pressure);
	strcat(tractionSlope,"/");
	strcat(tractionSlope,length);
	
	// viscosity
	int vexp = pexp+timeExp;
	if(unitsType==USER_UNITS)
		strcpy(viscosity,"F-T/L^2");
	else
	{	GetPrefix(vexp,"Pa-s",viscosity);
		if(strlen(viscosity)==0)
		{	GetPrefix(vexp+1,"P",viscosity);
			if(strlen(viscosity)==0)
				sprintf(viscosity,"10^%dPa-s",vexp);
		}
	}
	
	// heat capacity
	strcpy(heatCapacity,energy);
	strcat(heatCapacity,"/(");
	strcat(heatCapacity,mass);
	strcat(heatCapacity,"-K)");
	
	// power
	int powerexp = enerexp-timeExp;
	if(unitsType==USER_UNITS)
		strcpy(power,"F-L/T");
	else
	{	GetPrefix(powerexp,"W",power);
		if(strlen(power)==0)
		{	GetPrefix(powerexp+7,"(erg/s)",power);
			if(strlen(power)==0)
				sprintf(power,"10^%dW",powerexp);
		}
	}
	
	// conductivity
	if(unitsType==USER_UNITS)
		strcpy(conductivity,"F-L/(m-s-K)");
	else
	{	strcpy(conductivity,power);
		strcat(conductivity,"/(");
		strcat(conductivity,length);
		strcat(conductivity,"-K)");
	}

	// convection
	if(unitsType==USER_UNITS)
		strcpy(convection,"F-L/(m^2-s-K)");
	else
	{	strcpy(convection,power);
		strcat(convection,"/(");
		strcat(convection,length);
		strcat(convection,"^2-K)");
	}
	
	// heat flux
	if(unitsType==USER_UNITS)
		strcpy(heatFlux,"F-L/(m^2-s)");
	else
	{	strcpy(heatFlux,power);
		strcat(heatFlux,"/");
		strcat(heatFlux,length);
		strcat(heatFlux,"^2");
	}
	
	// bcarg
	strcpy(bcarg,timeu);
	strcat(bcarg,"/");
	strcat(bcarg,timeu);
	strcat(bcarg,"^-1");
	
	// concentration flux
	strcpy(concFlux,mass);
	strcat(concFlux,"/(");
	strcat(concFlux,length);
	strcat(concFlux,"^-2-");
	strcat(concFlux,timeu);
	strcat(concFlux,")");
	
	return true;
}

// get power of 3 label
void UnitsController::GetPrefix(int uexp,const char *SIUnit,char *makeUnits)
{
	if(uexp==15)
		strcpy(makeUnits,"P");
	else if(uexp==12)
		strcpy(makeUnits,"T");
	else if(uexp==9)
		strcpy(makeUnits,"G");
	else if(uexp==6)
		strcpy(makeUnits,"M");
	else if(uexp==3)
		strcpy(makeUnits,"k");
	else if(uexp==0)
		strcpy(makeUnits,"");
	else if(uexp==-2)
		strcpy(makeUnits,"c");
	else if(uexp==-3)
		strcpy(makeUnits,"m");
	else if(uexp==-6)
		strcpy(makeUnits,"u");
	else if(uexp==-9)
		strcpy(makeUnits,"n");
	else if(uexp==-12)
		strcpy(makeUnits,"p");
	else if(uexp==-15)
		strcpy(makeUnits,"f");
	else if(uexp==-15)
		strcpy(makeUnits,"a");
	else
	{	strcpy(makeUnits,"");
		return;
	}
	strcat(makeUnits,SIUnit);
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


