/********************************************************************************
    RigidMaterial.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Nov 15 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef RIGIDBCMATERIAL

#define RIGIDBCMATERIAL 11
#define RIGIDCONTACTMATERIAL 35

// same as constants without "CONTROL_" defined in boundary conditinos
#define CONTROL_X_DIRECTION 1
#define CONTROL_Y_DIRECTION 2
#define CONTROL_Z_DIRECTION 4
#define CONTROL_ANY_DIRECTION 7

#define CONTROL_TEMPERATURE 16
#define CONTROL_CONCENTRATION 32

// means an actual material, but only possible in multimaterial mode
#define RIGID_MULTIMATERIAL_MODE 8

class Expression;

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
        RigidMaterial(char *,int,int);
        
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
		void SetSettingFunction(char *,int);
        void ClearFunctions(void);
		void SetControlVelocity(double,int);
		void SetControlVelocity(Vector * Velocity);
	
		// Rigid methods subclassed might override
		virtual bool IsConstantVelocity(void) const;
		virtual bool GetVectorSetting(Vector *,bool *,double,Vector *) const;

		// accessors
		virtual const char *MaterialType(void) const;
        virtual double WaveSpeed(bool,MPMBase *) const;
		virtual bool IsRigid(void) const;			// override base class to return true
		virtual bool IsRigidBC(void) const;		// override base class to return true if appropriate
		virtual bool IsRigidContact(void) const;	// override base class to return true if appropriate
		virtual bool IsRigidBlock(void) const;		// override base class to return true if appropriate
	
	protected:
		Expression *function;
		Expression *function2;
		Expression *function3;
		Expression *Vfunction;

		bool useControlVelocity;
		int controlDirection;
		double controlVelocity;
		Vector ControlVelocityVector;

		static double varTime,xPos,yPos,zPos,delTime;
};

#endif
