/********************************************************************************
    JohnsonCook.hpp
    NairnMPM
    
    Created by John Nairn, August 12, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		Isoplasticity.hpp (IsotropicMat.hpp, Elastic.hpp, MaterialBase.hp)
********************************************************************************/

#ifndef JOHNSONCOOK

#define JOHNSONCOOK 16

#define YT_HISTORY 1
#define EPDOT_HISTORY 2

#include "Materials/IsoPlasticity.hpp"

class JohnsonCook : public IsoPlasticity
{
    public:

        // constructors and destructors
        JohnsonCook();
        JohnsonCook(char *);
        
        // initialize
		virtual char *InputMat(char *,int &);
		virtual void PrintMechanicalProperties(void);
		virtual void PrintYieldProperties(void);
		virtual void InitialLoadMechProps(int,int);
		virtual char *MaterialData(void);
				
		// methods
		virtual void LoadMechanicalProps(MPMBase *,int);
				
		// override plastic potential functions
		virtual double GetYield(MPMBase *,int,double);
 		virtual double GetKPrime(MPMBase *,int,double);
		virtual double GetK2Prime(MPMBase *,double,double);
		virtual void UpdatePlasticInternal(MPMBase *,int);
		
		virtual double GetHistory(int,char *);
       
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
		
    protected:
		double Bjc,Cjc,njc,ep0jc,Tmjc,mjc;
		double Bred,TjcTerm,edotMin,eminTerm;

};

#endif
