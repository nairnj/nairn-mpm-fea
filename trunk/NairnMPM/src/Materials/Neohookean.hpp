/********************************************************************************
	Neohookean.hpp
	nairn-mpm-fea

	Created by John Nairn on Sun Mar 10 2014.
	Copyright (c) 2014 John A. Nairn, All rights reserved.

	Neohookean hyperelastic material as an alternate to the Nehookean
		option in Mooney-Rivlin
 
	Dependencies
		HyperElastic.hpp
 ********************************************************************************/

#ifndef NEOHOOKEAN

#define NEOHOOKEAN 28

#include "Materials/HyperElastic.hpp"

// neohookean properties
typedef struct {
	double Gsp;
	double Lamesp;
	double Ksp;
} NeohookeanProperties;

class Neohookean : public HyperElastic
{
	public:
	
		// constructors and destructors
		Neohookean();
		Neohookean(char *);
	
		// initialize
		virtual char *InputMaterialProperty(char *,int &,double &);
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
		double G;
		double Etens,nu;
		double Lame;
		// double aI,betaI		// isotropic expanion defined in super classes
	
		NeohookeanProperties pr;
		double gamma0;
		int J_History;
};

#endif
