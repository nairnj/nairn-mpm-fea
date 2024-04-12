/********************************************************************************
    TIViscoelastic.hpp
    nairn-mpm-fea
    
    Created by John Nairn, Jan 7, 2021.
    Copyright (c) 2021 John A. Nairn, All rights reserved.

	Dependencies
        TransIsotropic.hpp (Elastic.hpp MaterialBase.hpp)
********************************************************************************/

#ifndef TIVISCOELASTIC1

#define TIVISCOELASTIC1 5
#define TIVISCOELASTIC2 6

#include "Materials/TransIsotropic.hpp"

enum { NO_PROP=0,GT_SERIES,GA_SERIES,KT_SERIES,N_SERIES,L_SERIES };

class TIViscoelastic : public TransIsotropic
{
    public:
        
        // constructors and destructors
        TIViscoelastic(char *,int);
		
		// initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
        void switchSeries(int,int);
        char *addToPkSeries(double **,int &,double &);
        char *addToTaukSeries(double **,int &);
        virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;
        void PrintPronySeries(const char *,double,double,double *,double *) const;
    
        // history
        char *InitHistoryData(char *,MPMBase *);
       virtual int NumberOfHistoryDoubles(void) const;

		// contitutive law methods
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int) const;
        void GetAlphaArgs(double,double,double &,double &) const;
        double GetPSArg(double,double) const;

		// accessors
  		virtual const char *MaterialType(void) const;
		virtual double WaveSpeed(bool,MPMBase *) const;
        virtual bool SupportsDiffusion(void) const;
		
    protected:
        double GT0,GA0,KT0,n0,ell0;
        int ntaus,whichOne;
        double *GTk,*tauGTk,GTe;
        double *GAk,*tauGAk,GAe;
        double *KTk,*tauKTk,KTe;
        double *nk,*taunk,ne;
        double *ellk,*tauellk,elle;
        int ntauGT,ntauGA,ntauKT,ntaun,ntauell;
        int currentPk,currentTauk;
        double betaAconc,betaTconc;
        int numHistory;
        double Tref,C1,C1base10,C2;
        double mref,Cm1,Cm2,Cm1base10,Cm2base10;
#ifdef OSPARTICULAS
        double tauMS,Cms;
#endif
};

#endif

