/********************************************************************************
    RigidMaterial.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Nov 15 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef RIGIDMATERIAL

#define RIGIDMATERIAL 11

// same as constants without "CONTROL_" defined in boundary conditinos
#define CONTROL_X_DIRECTION 1
#define CONTROL_Y_DIRECTION 2
#define CONTROL_Z_DIRECTION 4
#define CONTROL_ANY_DIRECTION 7

#define CONTROL_TEMPERATURE 16
#define CONTROL_CONCENTRATION 32

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
		static bool someSetTemperature;
		int mirrored;
	
        // constructors and destructors
        RigidMaterial();
        RigidMaterial(char *);
        
        // initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
        virtual void PrintMechanicalProperties(void) const;
		virtual int SetField(int,bool,int,int &);
		
		// RigidMaterial only methods
		bool RigidDirection(int) const;
		bool RigidTemperature(void) const;
		bool RigidConcentration(void) const;
		int SetDirection(void) const;
		bool GetValueSetting(double *,double,Vector *) const;
		bool GetVectorSetting(Vector *,bool *,double,Vector *) const;
		void SetSettingFunction(char *,int);
		void ReplaceSettingFunction(char *,int);
		bool IsConstantVelocity(void);
		void SetControlVelocity(double,int);
	
		// accessors
		virtual const char *MaterialType(void) const;
        virtual double WaveSpeed(bool,MPMBase *) const;
		virtual bool Rigid(void) const;			// override base class to return true
		virtual short RigidBC(void) const;		// override base class to return true if appropriate
		virtual short RigidContact(void) const;	// override base class to return true if appropriate
		
	protected:
		ROperation *function;
		ROperation *function2;
		ROperation *function3;
        ROperation *Vfunction;
		bool useControlVelocity;
		int controlDirection;
		double controlVelocity;
		static double varTime,xPos,yPos,zPos,delTime;
};

#endif
