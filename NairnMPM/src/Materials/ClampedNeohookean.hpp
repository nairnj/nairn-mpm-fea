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

#define JP_HISTORY 1
#define ELASTIC_DISNEY 0
#define ELASTIC_NH 1

class ClampedNeohookean : public Neohookean
{
	public:
		double critComp,critTens,hardening;
	int elasticModel;
	
		// constructors and destructors
		ClampedNeohookean();
		ClampedNeohookean(char *);
	
		// initialize
		virtual char *InputMat(char *,int &);
		virtual void PrintMechanicalProperties(void) const;
		virtual const char *VerifyAndLoadProperties(int);
		virtual void ValidateForUse(int) const;
		virtual char *InitHistoryData(void);
	
		// methods
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
	
		// accessors
		virtual const char *MaterialType(void) const;
		virtual int MaterialTag() const;
		virtual double GetHistory(int,char *) const;
		virtual bool SupportsArtificialViscosity(void) const;
		virtual double CurrentWaveSpeed(bool,MPMBase *) const;
	
	protected:
		double lamMin2,lamMax2;
	
};

#endif
