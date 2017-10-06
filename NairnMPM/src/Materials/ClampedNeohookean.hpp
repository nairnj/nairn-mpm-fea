/********************************************************************************
	ClampedNeohookean.hpp
	nairn-mpm-fea
	 
	Created by John Nairn on 2/3/15.
	Copyright (c) 2015 John A. Nairn, All rights reserved.
	 
	Neohookean hyperelastic material with plasticity based
	on maxima in tensile and compressive elongation
	Analogous t snow model used in Disney MPM paper
	 
	Dependencies
		Neohookean.hpp, HyperElastic.hpp
********************************************************************************/

#ifndef CLAMPEDNEOHOOKEAN

#define CLAMPEDNEOHOOKEAN 29

#include "Materials/Neohookean.hpp"

#define JP_HISTORY 2
#define ELASTIC_DISNEY 0
#define ELASTIC_NH 1

class ClampedNeohookean : public Neohookean
{
	public:
	
		// constructors and destructors
		ClampedNeohookean();
		ClampedNeohookean(char *);
	
		// initialize
		virtual char *InputMaterialProperty(char *,int &,double &);
		virtual void PrintMechanicalProperties(void) const;
		virtual const char *VerifyAndLoadProperties(int);
		virtual void ValidateForUse(int) const;
	
		// history data
		virtual int SizeOfHistoryData(void) const;
		virtual char *InitHistoryData(char *,MPMBase *);
		virtual double GetHistory(int,char *) const;
	
		// methods
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int) const;
	
		// accessors
		virtual const char *MaterialType(void) const;
		virtual bool SupportsArtificialViscosity(void) const;
		virtual double CurrentWaveSpeed(bool,MPMBase *,int) const;
	
	protected:
		double critComp,critTens,hardening;
		int elasticModel;
    	bool omitClamping;
	
		double lamMin2,lamMax2;
	
};

#endif
