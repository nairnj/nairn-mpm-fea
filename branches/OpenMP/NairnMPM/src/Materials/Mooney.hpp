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
        double G1,G2;
        bool rubber;
		// double aI,betaI		// isotropic expanion defined in super classes
 
        // constructors and destructors
        Mooney();
        Mooney(char *);
        
        // initialize
        virtual char *InputMat(char *,int &);
		virtual const char *VerifyProperties(int);
        virtual void InitialLoadMechProps(int,int);
		virtual void PrintMechanicalProperties(void);
        virtual char *InitHistoryData(void);
 		
		// methods
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int);
    
		// accessors
        virtual Tensor GetStress(Tensor *,double);
		virtual double WaveSpeed(bool,MPMBase *);
		virtual double ShearWaveSpeed(bool,MPMBase *);
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
        virtual double GetHistory(int,char *);
        virtual bool SupportsArtificialViscosity(void);
		
    protected:
		double G1sp, G2sp;
};

#endif
