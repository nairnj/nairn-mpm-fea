/********************************************************************************
	UnitsController.hpp
	NairnFEA/NairnMPM

	Created by John Nairn on 3/15/15.
	Copyright (c) 2015 John A. Nairn, All rights reserved.
********************************************************************************/

#ifndef _UNITSCONTROLLER_
#define _UNITSCONTROLLER_

enum { LEGACY_UNITS=0,CONSISTENT_UNITS,USER_UNITS };

enum { ERR_UNITS=0,CULENGTH_UNITS,TRACTIONSLOPE_UNITS,PRESSURE_UNITS,HEATCAPACITY_UNITS,
		C0_UNITS,VISCOSITY_UNITS,CONDUCTIVITY_UNITS,DENSITY_UNITS,STRESSINTENSITY_UNITS,
		FEAFORCE_UNITS,FEAWORK_UNITS,INTERFACEPARAM_UNITS,TIME_UNITS,CUMASS_UNITS,TARGETKE_UNITS,
		CUVELOCITY_UNITS,BCARG_UNITS,ALTTIME_UNITS,BCHEATFLUX_UNITS,BCCONCFLUX_UNITS,
		ALTVELOCITY_UNITS,CONVECTION_UNITS,HEATFUSION_UNITS,DIFFUSION_UNITS };

class UnitsController
{
	public:

		//  Constructors and Destructor
		UnitsController();
	
		// class methods
		static void OutputUnits(void);
		static char *ScaledPtr(char *,double &,double);
		static double Scaling(double);
		static const char *Label(int);
		static double UnitsAttribute(char *,int);
		static bool SetConsistentUnits(char *,char *,char *);
		static void GetPrefix(int,const char *,char *);
	
	private:
		static int unitsType;
	
		static char length[3],mass[3],timeu[3];
		static int lengthExp,massExp,timeExp;
		static char velocity[10],force[10],pressure[10],energy[10],density[10];
		static char errRate[15],tractionSlope[15],viscosity[15],heatFusion[18];
		static char heatCapacity[18],power[15],conductivity[20],convection[20];
		static char heatFlux[15],stressIntensity[20],bcarg[10],concFlux[20];
		static char diffusionTens[20];
};

#endif
