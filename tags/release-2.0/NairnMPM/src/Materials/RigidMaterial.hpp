/********************************************************************************
    RigidMaterial.hpp
    NairnMPM
    
    Created by John Nairn on Nov 15 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef RIGIDMATERIAL

#define RIGIDMATERIAL 11

// same aas constants without "CONTROL_" defined in boundary conditinos
#define CONTROL_X_DIRECTION 1
#define CONTROL_Y_DIRECTION 2
#define CONTROL_Z_DIRECTION 4

// means an actual material, but only possible in multimaterial mode
#define RIGID_MULTIMATERIAL_MODE 8

class ROperation;

#include "Materials/MaterialBase.hpp"

class RigidMaterial : public MaterialBase
{
    public:
		int setDirection;
		bool setTemperature;
		bool setConcentration;
		
        // constructors and destructors
        RigidMaterial();
        RigidMaterial(char *);
        
        // initialize
        virtual char *InputMat(char *,int &);
        virtual void PrintMechanicalProperties(void);
		virtual int SetField(int,bool,int);
		virtual void PreliminaryMatCalcs(void);
		
		// define abstract methods even though not used
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,int);
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int);
		
		// RigidMaterial only methods
		bool RigidDirection(int);
		bool RigidTemperature(void);
		bool RigidConcentration(void);
		bool GetSetting(double *,double);
		void SetSettingFunction(char *);
		
		// accessors
		virtual const char *MaterialType(void);
        virtual double WaveSpeed(bool);
		virtual short Rigid(void);			// override base class to return true
		virtual short RigidBC(void);		// override base class to return true if appropriate
		virtual int MaterialTag();
		
	protected:
		ROperation *function;
		static double varTime;
};

#endif