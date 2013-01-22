/********************************************************************************
    Viscoelastic.hpp
    NairnMPM
    
    Created by John Nairn on Wed Feb 5 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef VISCOELASTIC

#define VISCOELASTIC 7

#include "Materials/MaterialBase.hpp"

enum {XX_HISTORY=0,YY_HISTORY,XY_HISTORY,ZZ_HISTORY,XZ_HISTORY,YZ_HISTORY};

enum {VK_PROP=0,VA_PROP,VISCO_PROPS};

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
		virtual const char *VerifyProperties(int);
		virtual void InitialLoadMechProps(int,int);
		virtual void PrintMechanicalProperties(void);
		virtual void ValidateForUse(int);
		virtual char *InitHistoryData(void);
    
		// methods
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,int);
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int);
		
		// accessors
        virtual double WaveSpeed(bool,MPMBase *);
		virtual const char *MaterialType();
		virtual double GetHistory(int ,char *);
		virtual int MaterialTag();
		
    private:
        int currentGk,currentTauk;
		double Ge,dGe,CTE,CME,Ke;
        char read[VISCO_PROPS];
        
};

#endif

