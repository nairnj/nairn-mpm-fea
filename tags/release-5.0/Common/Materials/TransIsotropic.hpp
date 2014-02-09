/********************************************************************************
    TransIsotropic.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Tues Jan 28 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.

	Dependencies
		Elastic.hpp
********************************************************************************/

#ifndef TRANSISO1

#define TRANSISO1 2
#define TRANSISO2 3

#include "Materials/Elastic.hpp"

enum {EA_PROP,ET_PROP,GA_PROP,GT_PROP,NUT_PROP,TRANS_PROPS};
        
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
        virtual const char *VerifyAndLoadProperties(int);
        virtual void PrintMechanicalProperties(void) const;
#ifdef MPM_CODE
        virtual void FillTransportProperties(TransportProperties *);
		virtual void PrintTransportProperties(void) const;
#endif
		
		// methods
		void FillElasticProperties2D(ElasticProperties *,int,double,int) const;
#ifdef MPM_CODE
		void FillElasticProperties3D(MPMBase *,ElasticProperties *,int) const;
        virtual int SizeOfMechanicalProperties(int &) const;
		virtual void *GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *,void *) const;
		virtual void GetTransportProps(MPMBase *,int,TransportProperties *) const;
#else
		virtual void LoadMechanicalPropertiesFEA(int,double,int);
#endif
       
	   // accessors
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag() const;
#ifdef MPM_CODE
        virtual double WaveSpeed(bool,MPMBase *) const;
        virtual double MaximumDiffusion(void) const;
        virtual double MaximumDiffusivity(void) const;
		virtual double GetDiffZ(void) const;
		virtual double GetKcondZ(void) const;
#endif
		
	protected:
#ifdef FEA_CODE
		double lastMatAngle;
		int hasMatProps;
#endif
	
	private:
		int tiType;
};

#endif
