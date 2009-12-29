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
        double G0,K,aI,betaI;
        int ntaus;
        double *Gk,*tauk;
        
        // constructors and destructors
        Viscoelastic();
        Viscoelastic(char *);
        
        // methods
        virtual char *InputMat(char *,int &);
		virtual const char *VerifyProperties(int);
		virtual void InitialLoadMechProps(int,int);
		virtual void PrintMechanicalProperties(void);
        virtual char *MaterialData(void);
		virtual void ValidateUse(int);
    
		// methods
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,int);
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int);
		
		// accessors
        virtual double WaveSpeed(void);
		virtual const char *MaterialType();
		virtual double GetHistory(int ,char *);
		virtual int MaterialTag();
		
    private:
        int currentGk,currentTauk;
		double Ge,dGe,CTE,CME,Ke;
        char read[VISCO_PROPS];
        
};

#endif

