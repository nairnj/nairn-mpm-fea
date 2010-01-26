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
		virtual int MaterialTag();
		
	protected:
		ROperation *function;
		static double varTime;
};

#endif
