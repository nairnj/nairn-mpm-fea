/********************************************************************************
    MicrostructureModel.hpp
    NairnMPM
    
    Created by Tim Fagan in May 2012.
    Copyright (c) 2005 John A. Nairn, All rights reserved.

	Dependencies
		Isoplasticity.hpp (IsotropicMat.hpp, Elastic.hpp, MaterialBase.hp)
********************************************************************************/

#ifndef MICROSTRUCTUREMODEL

#define MICROSTRUCTUREMODEL 29
#define SQRT_ONETHIRD 0.5773502692
#define SQRT_THREE 1.732050808

#define YT_HISTORY 1
#define RHOC 2
#define RHOW 3
#define DSIZE 4
#define TDL 5
//#define EPDOT_HISTORY 6
//#define FR 7

#include "IsoPlasticity.hpp"

class MicrostructureModel : public IsoPlasticity
{
    public:
        // Constants
		//dislocation density in cell walls, dislocation density in cell, etc.
        double rhoW,rhoC,fLim,fo,fsto,SHM,N,alp,burg,K1,esal,esbe,disk1,tayM,sto;
		double fr0,tdl0,dSize0,rhoC0,rhoW0, MMG;
		
		
		
        // constructors and destructors
        MicrostructureModel();
        MicrostructureModel(char *);
        
        // initialize
        virtual char *InputMat(char *,int &);
		virtual void InitialLoadMechProps(int,int);
		virtual void PrintMechanicalProperties(void);
		virtual void PrintYieldProperties(void);
		virtual char *MaterialData(void);
				
		// override plastic potential functions
		virtual double GetYield(MPMBase *,int,double);
 		virtual double GetKPrime(MPMBase *,int,double);
		virtual double GetK2Prime(MPMBase *,double,double);
		virtual void UpdatePlasticInternal(MPMBase *,int);
       
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
		virtual double GetHistory(int,char *);
		
    protected:
        // specific values of yield strength & plastic modulus 
        double eqss;	// equivalent shear strain
		double rss;		// resolved shear strain
		double eqssra;	// equivalent shear streain rate;
		double rssra;	// resolved shear streain rate;
		double fr;		// volume fraction of cell walls;
		double tdl;		// total dislocation density;
		double dSize;	// cell size; <- save for later use?
		double cAdd,cRem,cDis;
		//double rhoCDot;
		double wAdd,wRem,wDis;
		//double rhoWDot;
		double rst;		// resolved shear stress
		double sigow;
		double rstw;	// shear stress in wall;
		double sigoc;
		double rstc;	// shear stress in cell;
		double yieldPrevious;	// previous yield stress
		double timeStep;	// the time increment for when archiving the strain rate.
		

};

#endif
