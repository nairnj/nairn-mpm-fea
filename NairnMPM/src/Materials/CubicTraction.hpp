/********************************************************************************
    CubicTraction.hpp
    nairn-mpm-fea
    
    Created by John Nairn on 4/6/08.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		TractionLaw.hpp, MaterialBase.hpp
********************************************************************************/

#ifndef CUBICTRACTIONMATERIAL

#define CUBICTRACTIONMATERIAL 14

#include "Materials/TractionLaw.hpp"

enum { CUBIC_DELN=0,CUBIC_DELT };

class CubicTraction : public TractionLaw
{
    public:
	
		// constructors and destructors
        CubicTraction(char *,int);
		
		// methods
		virtual const char *VerifyAndLoadProperties(int);
	
		// history data
        virtual char *InitHistoryData(char *,MPMBase *);
	
		// const methods
		virtual void PrintMechanicalProperties(void) const;
	
		// the traction law
		virtual void CrackTractionLaw(CrackSegment *,double,double,Vector *,Vector *,double);
		virtual double CrackWorkEnergy(CrackSegment *,double,double);
		virtual void CrackDissipatedEnergy(CrackSegment *,double &,double &);
		
		// accessors
		virtual const char *MaterialType(void) const;
		
	protected:
		double kI1,kII1;
};

#endif



