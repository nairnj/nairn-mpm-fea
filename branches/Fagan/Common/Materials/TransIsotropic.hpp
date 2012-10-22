/********************************************************************************
    TransIsotropic.hpp
    NairnMPM
    
    Created by John Nairn on Tues Jan 28 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.

	Dependencies
		Elastic.hpp
********************************************************************************/

#ifndef TRANSISO1

#define TRANSISO1 2
#define TRANSISO2 3

#include "Materials/Elastic.hpp"

enum {EA_PROP,ET_PROP,NUA_PROP,NUT_PROP,GA_PROP,
        AA_PROP,AT_PROP,GT_PROP,TRANS_PROPS};
        
enum {EX_PROP=0,EY_PROP,EZ_PROP,NUXY_PROP,NUYX_PROP,NUXZ_PROP,
        NUZX_PROP,NUYZ_PROP,NUZY_PROP,GXY_PROP,GXZ_PROP,
        GYZ_PROP,AX_PROP,AY_PROP,AZ_PROP,ORTHO_PROPS };

class TransIsotropic : public Elastic
{
    public:
        double EA,ET,nuA,nuT,GA,aA,aT,GT,KT,betaA,betaT;
        char read[ORTHO_PROPS];
#ifdef MPM_CODE
		double diffA,diffT,kcondA,kcondT;
#endif
        
        // constructors and destructors
        TransIsotropic();
        TransIsotropic(char *,int);
		
		// initialize
        virtual char *InputMat(char *,int &);
        virtual const char *VerifyProperties(int);
        virtual void PrintMechanicalProperties(void);
#ifdef MPM_CODE
		virtual void PrintTransportProperties(void);
#endif
		
		// methods
		virtual void LoadMechProps(int,double,int);
#ifdef MPM_CODE
		virtual void LoadMechanicalProps(MPMBase *,int);
		virtual void LoadTransportProps(MPMBase *,int);
#endif
       
	   // accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
#ifdef MPM_CODE
        virtual double WaveSpeed(bool,MPMBase *);
        virtual double MaximumDiffusion(void);
        virtual double MaximumDiffusivity(void);
		virtual double GetDiffZ(void);
		virtual double GetKcondZ(void);
#endif
		
	protected:
#ifdef MPM_CODE
		double lastTransAngle;
#endif
	
	private:
		int tiType;
};

#endif
