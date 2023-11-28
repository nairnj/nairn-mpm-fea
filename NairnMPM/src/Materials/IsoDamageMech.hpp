/********************************************************************************
	IsoDamageMech.hpp
	nairn-mpm-fea

	Created by John Nairn on Oct 6, 2021
	Copyright (c) 2021 John A. Nairn, All rights reserved.

	Dependencies
		IsoSoftening.hpp (IsotropicMat.hpp Elastic.hpp MaterialBase.hpp)
********************************************************************************/

#ifndef ISODAMAGEMECHANICS

#define ISODAMAGEMECHANICS 58

enum { STRESS_ENERGY_METRIC=0,TENSILE_ENERGY_METRIC,MIXED_MODE_METRIC,LAST_METRIC };

#include "Materials/IsoSoftening.hpp"

extern double predamageState,damageState;

class IsoDamageMech : public IsoSoftening
{
	public:
	
		// constructors and destructors
		IsoDamageMech(char *,int);
		
		// initialize
		virtual char *InputMaterialProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;
        virtual void SetInitialConditions(InitialCondition *,int,bool);

		// methods
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,
										int,Tensor *) const;
		virtual Tensor GetTrialStressIncrement(Tensor *,double,int) const;
		virtual Tensor GetDamageStressIncrement(Tensor *,double,int) const;
		virtual int ShouldInitiateFailure(Tensor *,double,Vector *,int) const;
		virtual double GetTractionFunction(Tensor *,Tensor *,int,double,double,double,double,double,double) const;
		virtual Tensor GetCe(Tensor *,int) const;

		// accessors
		virtual const char *MaterialType() const;

	protected:
		double Kred,Gred,Ered;
		int evolutionStyle;

};

#endif
