/********************************************************************************
    DDBHardening.hpp
	Written by Tim Fagan, June, 2013.
	
    nairn-mpm-fea
    
    Created by John Nairn, 11/12/2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.
 
	A dislocation density based polycrystal plasticity model
	See Toth et. al., J. Eng. Mat. and Techn., v124, 71-77 (2002)
	"Strain hardening at large strains as predicted by 
	dislocation based polycrystal plasticity model"

	Dependencies
		MaterialBase.hpp
********************************************************************************/

#ifndef _DDBHARDENING_

#define _DDBHARDENING_

#define SQRT_ONETHIRD 0.5773502692
#define SQRT_THREE 1.732050808

#define EP_HISTORY 0
#define YLDC_HISTORY 1 //5
#define EPDOT_HISTORY 2 //1
#define GRAIN_HISTORY 3 //4
#define RHOC_HISTORY 4 //2
#define RHOW_HISTORY 5 //3
#define RHOCDOT_HISTORY 6
#define RHOWDOT_HISTORY 7
#define NUMBER_HISTORY 8

#include "Materials/HardeningLawBase.hpp"

// plastic law properties
typedef struct {
	double SHM;
	double N;
	double fr;
	double rhoC;
	double rhoW;
	double rhoCtemp;
	double rhoWtemp;
	double rhoCDot;
	double rhoWDot;
	double rhoCDotTemp;
	double rhoWDotTemp;
	double rhoT;
	double dSize;
	double yieldP;
	double yieldInc;
	double yieldC;
	//double EVMYStress;
	//bool propertiesLoaded;
} DDBHProperties;


class DDBHardening : public HardeningLawBase
{
    public:
		
        
        // constructors and destructors
		DDBHardening();
		DDBHardening(MaterialBase *);
		
		// initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintYieldProperties(void) const;
		virtual void InitPlasticHistoryData(double *) const;
		
		
		// methods
		virtual int SizeOfHardeningProps(void) const;
        virtual void *GetCopyOfHardeningProps(MPMBase *,int,void *);
		virtual void DeleteCopyOfHardeningProps(void *,int) const;
		
		
		// hardening law core methods
        virtual double GetYield(MPMBase *,int,double,HardeningAlpha *,void *) const;
        virtual double GetKPrime(MPMBase *,int,double,HardeningAlpha *,void *) const;
        virtual double GetK2Prime(MPMBase *,double,double,HardeningAlpha *,void *) const;
		//virtual double GetYieldIncrement(MPMBase *,int,double,HardeningAlpha *,void *) const;
		//virtual double GetShearRatio(MPMBase *,double,double,void *) const;

		// Return Mapping
		virtual double SolveForLambdaBracketed(MPMBase *,int,double,Tensor *,double,double,double,double,HardeningAlpha *,void *) const;
		
		// History-dependent properties
		virtual int HistoryDoublesNeeded(void) const; 
		virtual void ElasticUpdateFinished(MPMBase *,int,double) const;
		virtual void UpdatePlasticInternal(MPMBase *,int,HardeningAlpha *) const;
		//virtual char *InitHistoryData(void);
		virtual double GetHistory(int,char *) const;
		
		
		// accessors
		virtual const char *GetHardeningLawName(void) const;


		
    protected:
	// properties independent of particle state
		double Atd, Btd;
		double fLim,fo,fsto,alp,burg,K1,esal,esbe,disk1,tayM,sto;
		double rhoC0,rhoW0, MMG, SHM0, N0;
		double dSize0,rhoT0;
		double tempDepend;
		// unique properties 
		/*double eqss;	// equivalent shear strain
		double rss;		// resolved shear strain
		double eqssra;	// equivalent shear streain rate;
		double rssra;	// resolved shear streain rate;
		double cAdd,cRem,cDis;
		double wAdd,wRem,wDis;
		double rst;		// resolved shear stress
		double sigow;
		double rstw;	// shear stress in wall;
		double sigoc;
		double rstc;	// shear stress in cell;
		*/

};

#endif

