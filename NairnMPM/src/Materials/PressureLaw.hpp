/********************************************************************************
	PressureLaw.hpp
	nairn-mpm-fea

	Created by John Nairn on 9/19/2013.
	Copyright (c) 2013 John A. Nairn, All rights reserved.

	Dependencies
		TractionLaw.hpp, MaterialBase.hpp
********************************************************************************/

#ifndef PRESSURELAWMATERIAL

#define PRESSURELAWMATERIAL 26

class Expression;

#include "Materials/TractionLaw.hpp"

class PressureLaw : public TractionLaw
{
	public:
	
		// constructors and destructors
		PressureLaw(char *,int);
	
		// methods
		virtual char *InputTractionLawProperty(char *,int &,double &);
		void SetStressFunction(char *);
		virtual void PrintMechanicalProperties(void) const;
		virtual void CalculateTimeFunction();
	
		// the traction law
		virtual void CrackTractionLaw(CrackSegment *,double,double,Vector *,Vector *,double);
		virtual double CrackTractionEnergy(CrackSegment *,double,double,bool);
	
		// accessors
		virtual const char *MaterialType(void) const;
	
	protected:
		Expression *function;
		double minCOD;
	
};

#endif
