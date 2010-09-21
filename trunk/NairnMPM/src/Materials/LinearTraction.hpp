/********************************************************************************
    LinearTraction.hpp
    NairnMPM
    
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
        LinearTraction(char *);
		
		// methods
		virtual const char *VerifyProperties(int);
		virtual void PrintMechanicalProperties(void);
	
		// the traction law
		virtual void CrackTractionLaw(CrackSegment *,double,double,double,double,double);
		virtual double CrackTractionEnergy(CrackSegment *,double,double,bool);
		
		// accessors
		virtual const char *MaterialType(void);
		virtual int MaterialTag();
};

#endif



