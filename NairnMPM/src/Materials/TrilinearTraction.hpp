/********************************************************************************
	TrilinearTraction.hpp
	nairn-mpm-fea

	Created by John Nairn on 9/17/10.
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
		CohesiveZone.hpp, TractionLaw.hpp, MaterialBase.hpp
 ********************************************************************************/

#ifndef TRILINEARTRACTIONMATERIAL

#define TRILINEARTRACTIONMATERIAL 20

#include "Materials/CohesiveZone.hpp"

class TrilinearTraction : public CohesiveZone
{
    public:
		
		// constructors and destructors
        TrilinearTraction(char *,int);
		
		// methods
		virtual char *InputTractionLawProperty(char *,int &,double &);
		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;
		
		// the traction law
		virtual void CrackTractionLaw(CrackSegment *,double,double,Vector *,Vector *,double);
		virtual double CrackTractionEnergy(CrackSegment *,double,double,bool);
		
		// accessors
		virtual const char *SetTLTractionLaw(double &,double &,double &,double &,double &,double &,double &);
		virtual const char *MaterialType(void) const;
		
	protected:
		double sI2,uI2;
		double sII2,uII2;
		bool break1is2I,break1is2II;
	
};

#endif



