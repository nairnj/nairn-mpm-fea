/********************************************************************************
    Mooney.hpp
    nairn-mpm-fea
    
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
		double Etens,nu;
        bool rubber;
		// double aI,betaI		// isotropic expanion defined in super classes
 
        // constructors and destructors
        Mooney();
        Mooney(char *);
        
        // initialize
        virtual char *InputMat(char *,int &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual char *InitHistoryData(void);
	
		// const methods
		virtual void PrintMechanicalProperties(void) const;
 		
		// methods
        virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
    
		// accessors
        virtual Vector ConvertJToK(Vector,Vector,Vector,int);
        virtual Tensor GetStress(Tensor *,double) const;
		virtual double WaveSpeed(bool,MPMBase *) const;
		virtual double ShearWaveSpeed(bool,MPMBase *) const;
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag() const;
        virtual double GetHistory(int,char *) const;
        virtual bool SupportsArtificialViscosity(void) const;
		virtual double GetCurrentRelativeVolume(MPMBase *) const;
		
    protected:
		double G1sp, G2sp;
		double gamma0;
		int J_History;
};

#endif
