/********************************************************************************
    LinearTraction.hpp
    nairn-mpm-fea
    
    Created by John Nairn on 4/1/08.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		TractionLaw.hpp, MaterialBase.hpp
********************************************************************************/

#ifndef LINEARTRACTIONMATERIAL

#define LINEARTRACTIONMATERIAL 13

#include "Materials/CohesiveZone.hpp"

class LinearTraction : public CohesiveZone
{
    public:
	
		// constructors and destructors
        LinearTraction(char *,int);
		
		// methods
		virtual const char *VerifyAndLoadProperties(int);
		virtual void PrintMechanicalProperties(void) const;
        virtual char *InitHistoryData(char *,MPMBase *);

		// the traction law
		virtual void CrackTractionLaw(CrackSegment *,double,double,Vector *,Vector *,double);
		virtual double CrackWorkEnergy(CrackSegment *,double,double);
		
		// accessors
		virtual const char *MaterialType(void) const;
};

#endif



