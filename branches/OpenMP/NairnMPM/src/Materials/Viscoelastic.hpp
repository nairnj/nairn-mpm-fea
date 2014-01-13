/********************************************************************************
    Viscoelastic.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef VISCOELASTIC

#define VISCOELASTIC 7

#include "Materials/MaterialBase.hpp"

enum {XX_HISTORY=0,YY_HISTORY,XY_HISTORY,ZZ_HISTORY,XZ_HISTORY,YZ_HISTORY};

class Viscoelastic : public MaterialBase
{
    public:
        double G0,K,aI;
        // double betaI;        // defined in superclass
        int ntaus;
        double *Gk,*tauk;
        
        // constructors and destructors
        Viscoelastic();
        Viscoelastic(char *);
        
        // innitialization methods
        virtual char *InputMat(char *,int &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual char *InitHistoryData(void);
	
		// const methods
		virtual void PrintMechanicalProperties(void) const;
		virtual void ValidateForUse(int) const;
    
		// methods
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,int,void *,ResidualStrains *) const;
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int,void *,ResidualStrains *) const;
        virtual double GetCpMinusCv(MPMBase *) const;
		
		// accessors
        virtual Vector ConvertJToK(Vector,Vector,Vector,int);
        virtual double WaveSpeed(bool,MPMBase *) const;
		virtual const char *MaterialType() const;
		virtual double GetHistory(int ,char *) const;
		virtual int MaterialTag() const;
		
    private:
        int currentGk,currentTauk;
		double Ge,dGe,CTE,CME,Ke,Ka2sp;
        
};

#endif

