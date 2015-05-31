/********************************************************************************
    TaitLiquid.hpp
    nairn-mpm-fea

    Created by John Nairn, Dec 4, 2013.
    Copyright (c) 2013 John A. Nairn, All rights reserved.

    Dependencies
        HyperElastic.hpp,MaterialBase.hpp
 ********************************************************************************/

#ifndef TAITLIQUID

#define TAITLIQUID 27
#define TAIT_C 0.0894

class ROperation;

#include "Materials/HyperElastic.hpp"

class TaitLiquid : public HyperElastic
{
    public:
    
        // constructors and destructors
        TaitLiquid();
        TaitLiquid(char *matName);
        
        // initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
        virtual void PrintMechanicalProperties(void) const;
        virtual void ValidateForUse(int) const;
        virtual void SetInitialParticleState(MPMBase *,int) const;
    
        // history variables
        char *InitHistoryData(void);
        double GetHistory(int num,char *historyPtr) const;
    
        // contitutive law methods
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
        
        // accessors
        virtual const char *MaterialType(void) const;
        virtual int MaterialTag() const;
        virtual double WaveSpeed(bool,MPMBase *) const;
        virtual double CurrentWaveSpeed(bool,MPMBase *) const;
        virtual Tensor GetStress(Tensor *sp,double pressure) const;
        virtual double GetCurrentRelativeVolume(MPMBase *) const;
		virtual void SetPressureFunction(char *);
    
    protected:
		// unique properties
		double viscosity;
	
        int J_history;
        double TwoEtasp;
		double gamma0;
		ROperation *function;
		static double xPos,yPos,zPos;
    
};

#endif

