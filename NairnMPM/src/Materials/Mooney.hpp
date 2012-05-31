/********************************************************************************
    Mooney.hpp
    NairnMPM
    
    Created by John Nairn on Fri Feb 28 2003.
    Copyright (c) 2003 John A. Nairn, All rights reserved.
 
	Moony-Rivlin hyperelastic material. A special case setting G2=0 is
	neo-Hookean material

	Dependencies
		HyperElastic.hpp
********************************************************************************/

#ifndef MOONEYRIVLIN

#define MOONEYRIVLIN 8

#include "Materials/HyperElastic.hpp"

enum {G1_PROP=0,G2_PROP,KBULK_PROP,CTE_PROP,MOONEY_PROPS};

class Mooney : public HyperElastic
{
    public:
        double G1,G2,Kbulk;
		// double aI,betaI		// isotropic expanion defined in super classes
 
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
        virtual double GetVolumetricTerms(double,double *);
    
		// accessors
		virtual double WaveSpeed(bool);
		virtual double ShearWaveSpeed(bool);
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
		
    private:
		double G1sp, G2sp, Ksp;
};

#endif
