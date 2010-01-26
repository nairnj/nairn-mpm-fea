/********************************************************************************
    Mooney.hpp
    NairnMPM
    
    Created by John Nairn on Fri Feb 28 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.

	Dependencies
		RubberElastic.hpp
********************************************************************************/

#ifndef MOONEYRIVLIN

#define MOONEYRIVLIN 8

#include "Materials/RubberElastic.hpp"

enum {C1_PROP=0,C2_PROP,CTE_PROP,MOONEY_PROPS};

class Mooney : public RubberElastic
{
    public:
        double C1,C2,aI,betaI;
 
        // constructors and destructors
        Mooney();
        Mooney(char *);
        
        // initialize
        virtual char *InputMat(char *,int &);
		virtual const char *VerifyProperties(int);
        virtual void InitialLoadMechProps(int,int);
		virtual void PrintMechanicalProperties(void);
 		
		// methods
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,int);
        virtual void MPMConstLaw(MPMBase *,double,double,double,double,double,double,double,double,double,double,int);
    
		// accessors
		virtual double WaveSpeed(bool);
		virtual bool ThreeDMaterial();
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
		
    private:
		double C1sp, C2sp;
        char read[MOONEY_PROPS];
};

#endif
